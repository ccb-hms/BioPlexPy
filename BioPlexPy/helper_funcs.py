import io
import requests
import anndata as ad
import pandas as pd
import networkx as nx

def getBioPlex(cell_line, version):
    '''
    Load BioPlex interactions data.
    
    This function loads BioPlex PPI data for
    for cell lines HEK293T and HCT116, note we
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

def bioplex2graph(bp_PPI_df):
    '''
    Convert BioPlex PPIs into a graph.
    
    This function converts representation of BioPlex PPIs into a graph data structure
    Representation of BioPlex PPIs in a NetworkX object from NetworkX. 

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame

    Returns
    -------
    NetworkX graph
        A NetworkX graph with Nodes = Uniprot Gene Symbols and Edges = interactions.

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> bp_293t_G = bioplex2graph(bp_293t_df) # (2) Turn the data into a graph
    '''
    # add isoform columns for Uniprot source & Uniprot target
    bp_PPI_df.loc[:,'isoformA'] = bp_PPI_df.UniprotA
    bp_PPI_df.loc[:,'isoformB'] = bp_PPI_df.UniprotB

    # reconstruct UniprotA/UniprotB columns without '-' isoform id
    UniprotA_new = []
    UniprotB_new = []
    for UniprotA, UniprotB in zip(bp_PPI_df.UniprotA, bp_PPI_df.UniprotB):

        if '-' in UniprotA:
            UniprotA_new.append(UniprotA.split('-')[0])
        else:
            UniprotA_new.append(UniprotA)

        if '-' in UniprotB:
            UniprotB_new.append(UniprotB.split('-')[0])
        else:
            UniprotB_new.append(UniprotB)

    # update columns for Uniprot source & Uniprot target to exclude isoform '-' ID
    bp_PPI_df.loc[:,'UniprotA'] = UniprotA_new
    bp_PPI_df.loc[:,'UniprotB'] = UniprotB_new
    
    # construct graph from BioPlex PPI data
    bp_G = nx.DiGraph()
    for source, target, pW, pNI, pInt in zip(bp_PPI_df.UniprotA, bp_PPI_df.UniprotB, bp_PPI_df.pW, bp_PPI_df.pNI, bp_PPI_df.pInt):
        bp_G.add_edge(source, target, pW=pW, pNI=pNI, pInt=pInt)
        
    # get mapping uniprot -> entrez & store as node attribute
    uniprot_entrez_dict = {}
    for uniprot_A, entrez_A in zip(bp_PPI_df.UniprotA, bp_PPI_df.GeneA):
        uniprot_entrez_dict[uniprot_A] = entrez_A
    for uniprot_B, entrez_B in zip(bp_PPI_df.UniprotB, bp_PPI_df.GeneB):
        uniprot_entrez_dict[uniprot_B] = entrez_B

    for node_i in bp_G.nodes():
        bp_G.nodes[node_i]["entrezid"] = uniprot_entrez_dict[node_i]

    # get mapping uniprot -> symbol & store as node attribute
    uniprot_symbol_dict = {}
    for uniprot_A, symbol_A in zip(bp_PPI_df.UniprotA, bp_PPI_df.SymbolA):
        uniprot_symbol_dict[uniprot_A] = symbol_A
    for uniprot_B, symbol_B in zip(bp_PPI_df.UniprotB, bp_PPI_df.SymbolB):
        uniprot_symbol_dict[uniprot_B] = symbol_B

    for node_i in bp_G.nodes():
        bp_G.nodes[node_i]["symbol"] = uniprot_symbol_dict[node_i]

    # get mapping uniprot -> isoform & store as node attribute
    uniprot_isoform_dict = {}
    for uniprot_A, isoform_A in zip(bp_PPI_df.UniprotA, bp_PPI_df.isoformA):
        uniprot_isoform_dict[uniprot_A] = isoform_A
    for uniprot_B, isoform_B in zip(bp_PPI_df.UniprotB, bp_PPI_df.isoformB):
        uniprot_isoform_dict[uniprot_B] = isoform_B

    for node_i in bp_G.nodes():
        bp_G.nodes[node_i]["isoform"] = uniprot_isoform_dict[node_i]
    
    return bp_G

def getCorum(complex_set = 'all', organism = 'Human'):
    '''
    Functionality for retrieving the CORUM protein complex data.

    Parameters
    ----------
    complex_set : str
        Takes input [all','core','splice'] (default 'all').
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