#!/usr/bin/env python

import pandas as pd
from pypdb import *
import itertools
from collections import Counter

def CORUM_to_PDB(Corum_DF, Complex_ID):
    '''
    Retreive PDB ID for protein structre corresponding to a CORUM complex.
    
    This function takes a CORUM complex ID and maps the corresponding 
    UniProt IDs for the complex to a PDB ID.

    Parameters
    ----------
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    PDB ID for CORUM complex

    Examples
    --------
    >>> Corum_DF = getCorum('core', 'Human') # (1) Obtain CORUM complexes
    >>> PDB_ID_Arp_2_3 = CORUM_to_PDB(Corum_DF, 27) # (2) Get PDB ID for specified protein complex (Arp 2/3 complex ID: 17)
    >>> PDB_ID_ING2 = CORUM_to_PDB(Corum_DF, 2851) # (3) Get PDB ID for specified protein complex (ING2 complex ID: 2851), demonstrates WARNINGS
    '''
    # get UniProt IDs for each protein in the CORUM complex
    uniprot_IDs_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(UniProt IDs)'].values[0].split(';')
    num_proteins_CORUM_complex_i = len(uniprot_IDs_complex_i)

    # Map from CORUM complex subunits given as UniProt IDs 
    # via [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html) to PDB structures:
    # "A summary of the UniProt to PDB mappings showing the UniProt accession followed by a semicolon-separated list of PDB four letter codes."
    uniprot_pdb_mapping_df = pd.read_csv("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_pdb.csv.gz", header = 1, sep = ',', compression = 'gzip')

    # set UniProt IDs as index
    uniprot_pdb_mapping_df.set_index('SP_PRIMARY', drop = True, inplace = True)

    # convert col PDB semicolon-separated list into Python list
    uniprot_pdb_mapping_df.loc[:,'PDB'] = [PDB_codes_i.split(';') for PDB_codes_i in uniprot_pdb_mapping_df.PDB]

    # get PDB IDs that map to each UniProt ID
    PDB_IDs_complex_i = []
    for uniprot_ID_i_complex_i in uniprot_IDs_complex_i:

        # check to see if UniProt ID exists in mapping
        if uniprot_ID_i_complex_i in uniprot_pdb_mapping_df.index:

            # append to list of PDB IDs, take ALL PDB IDs in mapped list
            mapped_PDB_ID_i = uniprot_pdb_mapping_df.loc[uniprot_ID_i_complex_i,:].values[0]
            PDB_IDs_complex_i.append(mapped_PDB_ID_i)

        else:
            print(f'WARNING: {uniprot_ID_i_complex_i} does not have any corresponding PDB IDs mapped.')

    # check to see there are any PDB IDs that are present in all Uniprot - PDB mappings
    PDB_IDs_complex_i = list(itertools.chain(*PDB_IDs_complex_i)) # flatten list of PDB IDs that mapped to UniProt IDs
    # convert to uppercase
    PDB_IDs_complex_i = [PDB_ID.upper() for PDB_ID in PDB_IDs_complex_i]
    PDB_IDs_complex_i_count = pd.Series(Counter(PDB_IDs_complex_i))

    # check if one PDB ID maps to all UniProt IDs
    if sum(PDB_IDs_complex_i_count == len(uniprot_IDs_complex_i)) == 1:

        # use this PDB ID for CORUM complex
        PDB_ID_structure_for_CORUM_complex_i = list(PDB_IDs_complex_i_count[PDB_IDs_complex_i_count == len(uniprot_IDs_complex_i)].index)[0]

    else:

        # check to see if multiple PDB IDs map to all UniProt IDs, then search through those
        if sum(PDB_IDs_complex_i_count == len(uniprot_IDs_complex_i)) > 1:

            # get PDB IDs that map to all Uniprot IDs
            PDB_IDs_mapped_to_complex_i = list(PDB_IDs_complex_i_count[PDB_IDs_complex_i_count == len(uniprot_IDs_complex_i)].index)

        # if no PDB ID maps to all UniProt IDs, 
        # then run a query with all unique PDB IDs that mapped to any UniProt ID and search the resulting PDB IDs to find complex 
        else:
            PDB_IDs_complex_i = list(set(PDB_IDs_complex_i)) # if redundant PDB IDs
            PDB_IDs_complex_i = ' '.join(PDB_IDs_complex_i) # convert seperate PDB ids into one string

            # use PDB IDs from SIFTS mapping as search terms to find associated complex PDB ID in PDB database
            PDB_IDs_mapped_to_complex_i = Query(PDB_IDs_complex_i).search()
            
        # if query returns an empty list, raise warning
        if len(PDB_IDs_mapped_to_complex_i) == 0:
            print(f'WARNING: Could not map PDB ID to this CORUM complex ID.')
            PDB_ID_structure_for_CORUM_complex_i = None

        else:
            # iterate through PDB ID & retreive metadata
            PDB_protein_count = [] # stores number of polymer proteins for this structure
            PDB_deposit_date = [] # store the date of deposit for this structure
            for PDB_ID in PDB_IDs_mapped_to_complex_i:

                # retreive metadata for this structure from PDB
                PDB_structure_all_info = get_info(PDB_ID)
                PDB_protein_count.append(PDB_structure_all_info['rcsb_entry_info']['polymer_entity_count_protein'])
                PDB_deposit_date.append(pd.to_datetime(PDB_structure_all_info['rcsb_accession_info']['deposit_date']))

            # convert CORUM complex i - associated PDB IDs into DataFrame w/ # proteins & resolution
            complex_i_PDBs_df = pd.DataFrame(index = PDB_IDs_mapped_to_complex_i)
            complex_i_PDBs_df.loc[:,'num_proteins'] = PDB_protein_count
            complex_i_PDBs_df.loc[:,'deposit_date'] = PDB_deposit_date

            # column for number of proteins in PDB structure different from num proteins listed in CORUM complex
            complex_i_PDBs_df.loc[:,'num_proteins_PDB_CORUM_diff'] = abs(complex_i_PDBs_df.num_proteins - num_proteins_CORUM_complex_i)

            # pick PDB structure that has same number of proteins/chains as CORUM complex (or matches closest), then rank by most recent deposit date
            complex_i_PDBs_df.sort_values(by = ['num_proteins_PDB_CORUM_diff','deposit_date'], ascending = [True, False], inplace = True)
            PDB_ID_structure_for_CORUM_complex_i = complex_i_PDBs_df.index[0] # take PDB ID corresponding to top row after ranking
            
            # print a warning if the number of proteins for PDB ID complex differs from number of subunits in CORUM complex
            num_proteins_diff = complex_i_PDBs_df.loc[PDB_ID_structure_for_CORUM_complex_i, 'num_proteins_PDB_CORUM_diff']
            if num_proteins_diff != 0:
                print(f'WARNING: The number of proteins in PDB ID {PDB_ID_structure_for_CORUM_complex_i} and number of subunits in CORUM complex {Complex_ID} differs by {num_proteins_diff}.')
        
    return PDB_ID_structure_for_CORUM_complex_i