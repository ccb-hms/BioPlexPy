#!/usr/bin/env python

import io
import itertools

import anndata as ad
import pandas as pd
import requests
from pypdb import *


def getBioPlex(cell_line, version):
    '''
    Load BioPlex interactions data.
    
    This function loads BioPlex PPI data for
    cell lines HEK293T and HCT116, note we
    only have version 1.0 for HCT116 cells.
    
    Parameters
    ----------
    cell_line : str
        Takes input ['293T','HCT116'].
    version : str
        Takes input ['3.0','1.0','2.0'].
    
    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a PPI interaction.
    
    Examples
    --------
    >>> bp_293t = getBioPlex('293T', '1.0')
    >>> bp_293t.head(1)
       GeneA   GeneB UniprotA UniprotB SymbolA SymbolB            pW       pNI      pInt
    0    100  728378   P00813   A5A3E0     ADA   POTEF  2.605947e-09  0.000333  0.999667
    '''
    if (f'{cell_line}.{version}' not in 
        ['293T.1.0','293T.2.0','293T.3.0','HCT116.1.0']):
        print('dataset not available for this Cell Line - Version')
        
    else:
        if f'{cell_line}.{version}' == '293T.1.0':
            file_ext = 'interactionList_v2'
        elif f'{cell_line}.{version}' == '293T.2.0':
            file_ext = 'interactionList_v4a'
        elif f'{cell_line}.{version}' == '293T.3.0':
            file_ext = '293T_Network_10K_Dec_2019'
        elif f'{cell_line}.{version}' == 'HCT116.1.0':
            file_ext = 'HCT116_Network_5.5K_Dec_2019'

        BioPlex_interactions_df = pd.read_csv(
                f"https://bioplex.hms.harvard.edu/data/BioPlex_{file_ext}.tsv", 
                sep = '\t')
        
    # if pulling 293T cell line version 1.0 or 2.0, change column names to 
    # standardize across datasets for input into other functions
    if (cell_line == '293T') and (version == '1.0'):
        BioPlex_interactions_df.rename({'Gene A':'GeneA','Gene B':'GeneB',
            'Uniprot A':'UniprotA','Uniprot B':'UniprotB',
            'Symbol A':'SymbolA','Symbol B':'SymbolB','p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)

    if (cell_line == '293T') and (version == '2.0'):
        BioPlex_interactions_df.rename({'p(Wrong)':'pW',
            'p(No Interaction)':'pNI','p(Interaction)':'pInt'}, 
            axis = 1, inplace = True)
    
    return BioPlex_interactions_df

def getGSE122425():
    '''
    Retrieve HEK293 RNAseq expression data.
    
    Returns
    -------
    adata : AnnData object
        SummarizedExperiment of HEK293 raw count with an 
        added layer storing rpkm.
    
    Examples
    --------
    >>> HEK293_adata = getGSE122425()
    >>> HEK293_adata
    AnnData object with n_obs × n_vars = 57905 × 6
        obs: 'SYMBOL', 'KO', 'GO', 'length'
        layers: 'rpkm'
    >>> print(HEK293_adata.obs_names[:10].tolist())
    ['ENSG00000223972', 'ENSG00000227232', 'ENSG00000243485', 'ENSG00000237613', 'ENSG00000268020', 'ENSG00000240361', 'ENSG00000186092', 'ENSG00000238009', 'ENSG00000239945', 'ENSG00000233750']
    '''
    # specify URL where data is stored
    baseURL = ('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/'
               'GSE122425/suppl/')
    filename = 'GSE122425_all.counts.293_vs_293NK.edgeR_all.xls.gz'
    outFilePath = filename[:-3]

    # stream the file as bytes into memory using 
    # io.bytesIO and decompress using pandas
    response = requests.get(baseURL + filename)
    content = response.content
    GSE122425_df = pd.read_csv(io.BytesIO(content), 
                               sep='\t', compression='gzip')

    # create annot for observations (rows)
    obs = GSE122425_df.loc[:,['gene_id','GeneSymbol','KO','GO','length']]
    obs.set_index('gene_id',inplace=True)
    obs.rename(mapper={'GeneSymbol':'SYMBOL'}, inplace=True, axis=1)

    # we have raw counts and rpkms here in one matrix
    # raw counts
    raw_X = (GSE122425_df.loc[:,
                              ['HEK293NK-SEQ1','HEK293NK-SEQ2',
                               'HEK293NK-SEQ3','HEK293-SEQ1','HEK293-SEQ2',
                               'HEK293-SEQ3']].values)
    # annot for variables (cols)
    raw_var = pd.DataFrame(index=['NK.1','NK.2','NK.3','WT.1','WT.2','WT.3'])
    
    # convert to AnnData object (default datatype is 'float32')
    adata = ad.AnnData(raw_X, obs=obs, var=raw_var, dtype='int32')

    # store rpkms as a layer
    rpkm_X = (GSE122425_df.loc[:,['HEK293NK-SEQ1_RPKM','HEK293NK-SEQ2_RPKM',
                                  'HEK293NK-SEQ3_RPKM','HEK293-SEQ1_RPKM',
                                  'HEK293-SEQ2_RPKM',
                                  'HEK293-SEQ3_RPKM']].values)
    adata.layers["rpkm"] = rpkm_X

    return adata

def getCorum(complex_set = 'all', organism = 'Human'):
    '''
    Functionality for retrieving the CORUM protein complex data.

    Parameters
    ----------
    complex_set : str
        Takes input ['all','drug','splice','partial'] (default 'all').
        Maps to CORUM files: 'all' -> corum_allComplexes.txt, 
        'drug' -> corum_drugs.txt, 
        'splice' -> corum_spliceComplexes.txt,
        'partial' -> corum_partialComplexes.txt.
    organism : str
        Takes input ['Bovine','Dog','Hamster','Human','MINK','Mammalia',
        'Mouse','Pig','Rabbit','Rat'] (default 'Human').

    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a CORUM complex.

    Examples
    --------
    >>> CORUM_df = getCorum()
    >>> CORUM_df = getCorum('splice', 'Human')
    >>> CORUM_df.size > 0
    True
    >>> 'complex_id' in CORUM_df.columns
    True
    >>> 'subunits_uniprot_id' in CORUM_df.columns
    True
    '''
    # Map complex_set input to CORUM filename
    complex_set_map = {
        'all': 'corum_allComplexes.txt',
        'partial': 'corum_partialComplexes.txt',
        'drug': 'corum_drugs.txt',
        'splice': 'corum_spliceComplexes.txt'
    }
    
    if complex_set not in complex_set_map:
        raise ValueError(f"complex_set must be one of {list(complex_set_map.keys())}, "
                        f"got '{complex_set}'")
    
    # specify URL
    baseURL = 'https://zenodo.org/records/17419058/files/'
    filename = complex_set_map[complex_set]

    # stream the file as bytes into memory using
    # io.BytesIO and read using pandas
    response = requests.get(baseURL + filename)
    response.raise_for_status()  # raise exception for bad status codes
    content = response.content
    CORUM_df = pd.read_csv(io.BytesIO(content), sep='\t')

    # filter to keep only CORUM sets for a specific organism
    CORUM_df = CORUM_df[CORUM_df.organism == organism]
    CORUM_df.reset_index(inplace = True, drop = True)
    
    return CORUM_df

def get_UniProts_from_CORUM(Corum_DF, Complex_ID):
    '''
    Retrieve set of UniProt IDs corresponding to a CORUM complex ID.
    
    This function takes a CORUM complex ID and CORUM complex DataFrame
    and returns the corresponding UniProt IDs.

    Parameters
    ----------
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    UniProt IDs
        A list of UniProt IDs for the CORUM complex specified.

    Examples
    --------
    # (1) Obtain CORUM complexes
    # (2) Get set of UniProt IDs for specified protein 
    #     complex (Arp 2/3 complex ID: 27)
    >>> Corum_DF = getCorum('all', 'Human')
    >>> UniProts_Arp_2_3 = get_UniProts_from_CORUM(Corum_DF, Complex_ID = 27)
    >>> len(UniProts_Arp_2_3) > 0
    True
    >>> isinstance(UniProts_Arp_2_3, list)
    True
    '''
    # get UniProt IDs for each protein in the CORUM complex
        # NOTE: CORUM 5.1 now has a uniprot mapping txt file, not sure if better than DF?
    uniprot_IDs_list = (Corum_DF[Corum_DF.complex_id == Complex_ID].loc[:,
                                'subunits_uniprot_id'].values[0].split(';'))
    return uniprot_IDs_list

def get_PDB_from_UniProts(uniprot_IDs_list):
    '''
    Retreive PDB IDs for protein structures corresponding to set of UniProt IDs.
    
    This function takes a list of UniProt IDs and maps the corresponding 
    UniProt IDs (from the UniProt IDs input or CORUM complex ID) to PDB IDs 
    using the SIFTS project. Some metadata for each PDB ID is pulled from PDB 
    and stored in a DataFrame that is returned.

    Parameters
    ----------
    UniProt IDs : list

    Returns
    -------
    PDB IDs and associated metadata
        Pandas DataFrame of PDB IDs that map to the UniProt IDs input, 
        or corresponding UniProt IDs from the CORUM complex specified.

    Examples
    --------

    # this test is not consistent > PDB_ID_Arp_2_3.UniProts_mapped_to_PDB[1]
    ['O15144', 'P61158', 'P61160', 'O15145', 'P59998']
    # (1) Get set of PDB IDs for list of UniProt IDs that correspond to Arp 2/3
    >>> PDB_ID_Arp_2_3 = get_PDB_from_UniProts(['Q92747','O15144','P61158','P61160','O15145','P59998','O15511'])
    >>> PDB_ID_Arp_2_3.size
    15
    >>> type(PDB_ID_Arp_2_3)
    <class 'pandas.core.frame.DataFrame'>
    '''
    # get number of proteins in query
    num_proteins = len(uniprot_IDs_list)

    # Map from CORUM complex subunits given as UniProt IDs 
    # via [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html) 
    # to PDB structures:
    # "A summary of the UniProt to PDB mappings showing the UniProt accession 
    #  followed by a semicolon-separated list of PDB four letter codes."
    uniprot_pdb_mapping_df = pd.read_csv(
            ('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/'
            'flatfiles/csv/uniprot_pdb.csv.gz'), header = 1, 
             sep = ',', compression = 'gzip')

    # set UniProt IDs as index
    uniprot_pdb_mapping_df.set_index('SP_PRIMARY', drop = True, inplace = True)

    # convert col PDB semicolon-separated list into Python list
    uniprot_pdb_mapping_df.loc[:,'PDB'] = ([PDB_codes_i.split(';') 
                                for PDB_codes_i in uniprot_pdb_mapping_df.PDB])

    # get PDB IDs that map to each UniProt ID
    PDB_IDs_for_uniprot_dict = {}
    for uniprot_ID_i in uniprot_IDs_list:

        # check to see if UniProt ID exists in mapping
        if uniprot_ID_i in uniprot_pdb_mapping_df.index:

            # append to list of PDB IDs, take ALL PDB IDs in mapped list
            mapped_PDB_ID_i = (uniprot_pdb_mapping_df.loc[uniprot_ID_i,
                                                          :].values[0])

            # convert to uppercase
            mapped_PDB_ID_i = [PDB_ID.upper() for PDB_ID in mapped_PDB_ID_i]

            PDB_IDs_for_uniprot_dict[uniprot_ID_i] = mapped_PDB_ID_i

        else:
            print(f'WARNING: {uniprot_ID_i} does not have any '
                  'corresponding PDB IDs mapped.')

    # create dictionary of PDB IDs and store list of 
    # Uniprot IDs each one is mapped to
    # flatten list of PDB IDs that mapped to UniProt IDs
    unique_PDB_IDs = list(set(
            itertools.chain(*list(PDB_IDs_for_uniprot_dict.values()))))
    uniprot_IDs_list_for_PDB_dict = {}
    for PDB_ID in unique_PDB_IDs:

        uniprot_to_PDB_i = []
        for uniprot_ID in PDB_IDs_for_uniprot_dict.keys():
            if PDB_ID in PDB_IDs_for_uniprot_dict[uniprot_ID]:
                uniprot_to_PDB_i.append(uniprot_ID)

        uniprot_IDs_list_for_PDB_dict[PDB_ID] = uniprot_to_PDB_i

    uniprot_IDs_list_for_PDB_series = pd.Series(uniprot_IDs_list_for_PDB_dict)

    # if no PDB IDs mapped to UniProt IDs (empty list), raise warning
    if len(uniprot_IDs_list_for_PDB_series) == 0:
        print('WARNING: Could not map PDB ID to this CORUM '
              'complex ID or UniProt IDs.')
        complex_i_PDBs_df = None

    else:
        # iterate through PDB ID & retreive metadata
        # stores number of polymer proteins for this structure
        PDB_protein_count = []
        # store the date of deposit for this structure
        PDB_deposit_date = []
        # store the title of the citation for this structure
        PDB_citation_title = []
        # store list of UniProt IDs that mapped to each PDB ID
        PDB_to_uniprot_map_list = []
        for PDB_ID in uniprot_IDs_list_for_PDB_series.index:

            # retreive metadata for this structure from PDB
            PDB_structure_all_info = get_info(PDB_ID)
            PDB_protein_count.append(
                (PDB_structure_all_info['rcsb_entry_info']
                ['polymer_entity_count_protein']))
            PDB_deposit_date.append(
                (pd.to_datetime(PDB_structure_all_info['rcsb_accession_info']
                ['deposit_date'])))
            PDB_citation_title.append(
                PDB_structure_all_info['rcsb_primary_citation']['title'])
            PDB_to_uniprot_map_list.append(
                uniprot_IDs_list_for_PDB_series[PDB_ID])

        # convert CORUM complex i - associated PDB IDs into 
        # DataFrame w/ # proteins & resolution
        UniProt_assoc_PDBs_df = (pd.DataFrame(
            index = uniprot_IDs_list_for_PDB_series.index))
        UniProt_assoc_PDBs_df.loc[:,'num_proteins'] = PDB_protein_count
        UniProt_assoc_PDBs_df.loc[:,'deposit_date'] = PDB_deposit_date
        UniProt_assoc_PDBs_df.loc[:,'citation_title'] = PDB_citation_title
        # the UniProt IDs that mapped to this PDB ID from SIFTS
        (UniProt_assoc_PDBs_df.loc[:,
            'UniProts_mapped_to_PDB']) = PDB_to_uniprot_map_list

        # column for number of proteins in PDB structure different from
        # num proteins listed in CORUM complex or input by user
        (UniProt_assoc_PDBs_df.loc[:,
            'num_proteins_diff_btwn_PDB_and_UniProts_input']) = abs(
                UniProt_assoc_PDBs_df.num_proteins - num_proteins)

        # pick PDB structure that has same number of proteins/chains as CORUM
        # complex (or matches closest), then rank by most recent deposit date
        # tyl: recent deposit date fluctuates?
        UniProt_assoc_PDBs_df.sort_values(
            by = (['num_proteins_diff_btwn_PDB_and_UniProts_input',
                   'deposit_date']), 
            ascending = [True, False], inplace = True)
        
    return UniProt_assoc_PDBs_df