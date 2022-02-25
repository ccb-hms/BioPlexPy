#!/usr/bin/env python

import io
import requests
import anndata as ad
import pandas as pd

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
    >>> bp_hct116 = getBioPlex('HCT116', '1.0')
    '''
    if f'{cell_line}.{version}' not in ['293T.1.0','293T.2.0','293T.3.0','HCT116.1.0']:
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

        BioPlex_interactions_df = pd.read_csv(f"https://bioplex.hms.harvard.edu/data/BioPlex_{file_ext}.tsv", sep = '\t')
        
    # if pulling 293T cell line version 1.0 or 2.0, change column names to standardize across datasets for input into other functions
    if (cell_line == '293T') and (version == '1.0'):
        BioPlex_interactions_df.rename({'Gene A':'GeneA','Gene B':'GeneB','Uniprot A':'UniprotA','Uniprot B':'UniprotB','Symbol A':'SymbolA','Symbol B':'SymbolB','p(Wrong)':'pW','p(No Interaction)':'pNI','p(Interaction)':'pInt'}, axis = 1, inplace = True)

    if (cell_line == '293T') and (version == '2.0'):
        BioPlex_interactions_df.rename({'p(Wrong)':'pW','p(No Interaction)':'pNI','p(Interaction)':'pInt'}, axis = 1, inplace = True)
    
    return BioPlex_interactions_df

def getGSE122425():
    '''
    Retrieve HEK293 RNAseq expression data.
    
    Returns
    -------
    adata : AnnData object
        SummarizedExperiment of HEK293 raw count with an added layer storing rpkm.
    
    Examples
    --------
    >>> HEK293_adata = getGSE122425()
    '''
    # specify URL where data is stored
    baseURL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122425/suppl/'
    filename = 'GSE122425_all.counts.293_vs_293NK.edgeR_all.xls.gz'
    outFilePath = filename[:-3]

    # stream the file as bytes into memory using io.bytesIO and decompress using pandas
    response = requests.get(baseURL + filename)
    content = response.content
    GSE122425_df = pd.read_csv(io.BytesIO(content), sep='\t', compression='gzip')

    # create annot for observations (rows)
    obs = GSE122425_df.loc[:,['gene_id','GeneSymbol','KO','GO','length']]
    obs.set_index('gene_id',inplace=True)
    obs.rename(mapper={'GeneSymbol':'SYMBOL'}, inplace=True, axis=1)

    # we have raw counts and rpkms here in one matrix
    # raw counts
    raw_X = GSE122425_df.loc[:,['HEK293NK-SEQ1','HEK293NK-SEQ2','HEK293NK-SEQ3','HEK293-SEQ1','HEK293-SEQ2','HEK293-SEQ3']].values
    raw_var = pd.DataFrame(index=['NK.1','NK.2','NK.3','WT.1','WT.2','WT.3']) # annot for variables (cols)
    
    # convert to AnnData object (default datatype is 'float32')
    adata = ad.AnnData(raw_X, obs=obs, var=raw_var, dtype='int32')

    # store rpkms as a layer
    rpkm_X = GSE122425_df.loc[:,['HEK293NK-SEQ1_RPKM','HEK293NK-SEQ2_RPKM','HEK293NK-SEQ3_RPKM','HEK293-SEQ1_RPKM','HEK293-SEQ2_RPKM','HEK293-SEQ3_RPKM']].values
    adata.layers["rpkm"] = rpkm_X

    return adata

def getCorum(complex_set = 'all', organism = 'Human'):
    '''
    Functionality for retrieving the CORUM protein complex data.

    Parameters
    ----------
    complex_set : str
        Takes input ['all','core','splice'] (default 'all').
    organism : str
        Takes input ['Bovine','Dog','Hamster','Human','MINK','Mammalia','Mouse','Pig','Rabbit','Rat'] (default 'Human').

    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a CORUM complex.

    Examples
    --------
    >>> CORUM_df = getCorum()
    >>> CORUM_df = getCorum('core', 'Human')
    '''
    # specify URL where data is stored
    baseURL = 'https://mips.helmholtz-muenchen.de/corum/download/'
    filename = f'{complex_set}Complexes.txt.zip'
    outFilePath = filename[:-4]

    # stream the file as bytes into memory using io.bytesIO and decompress using pandas
    response = requests.get(baseURL + filename)
    content = response.content
    CORUM_df = pd.read_csv(io.BytesIO(content), sep='\t', compression='zip')

    # filter to keep only CORUM sets for a specific organism
    CORUM_df = CORUM_df[CORUM_df.Organism == organism]
    CORUM_df.reset_index(inplace = True, drop = True)
    
    return CORUM_df