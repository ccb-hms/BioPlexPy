import io
import requests
import anndata as ad
import pandas as pd

def getBioPlex(cell_line, version):

    '''
    Function to Load BioPlex interactions data
    
    cell_line: takes input ['293T','HCT116']
    version: takes input ['3.0','1.0','2.0']
    
    data for cell lines: HEK293T and HCT116
    we only have version 1.0 for HCT116 cells
    
    Column Descriptions
    GeneA: Entrez Gene ID for the first interacting protein
    GeneB: Entrez Gene ID for the second interacting protein
    UniprotA: Uniprot ID for the first interacting protein
    UniprotB: Uniprot ID for the second interacting protein
    SymbolA: Symbol for the first interacting protein
    SymbolB: Symbol for the second interacting protein
    p(Wrong ID): Probability of wrong protein ID (CompPASS-Plus)
    p(NotInteractor): Probability of nonspecific background (CompPASS-Plus)
    p(Interactor): Probability of high-confidence interaction (CompPASS-Plus)
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
    Function to retrieve HEK293 RNAseq expression data
    
    output is two AnnData objects: (1) HEK293 raw counts, (2) HEK293 rpkm
    example usage:
    
        HEK293_adata_raw, HEK293_adata_rpkm = getGSE122425()
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

    # we have raw counts and rpkms here in one matrix...
    # let's pull them out seperately and make each one an assay

    # raw counts
    raw_X = GSE122425_df.loc[:,['HEK293NK-SEQ1','HEK293NK-SEQ2','HEK293NK-SEQ3','HEK293-SEQ1','HEK293-SEQ2','HEK293-SEQ3']].values
    raw_var = pd.DataFrame(index=['NK.1','NK.2','NK.3','WT.1','WT.2','WT.3']) # annot for variables (cols)

    # rpkms
    rpkm_X = GSE122425_df.loc[:,['HEK293NK-SEQ1_RPKM','HEK293NK-SEQ2_RPKM','HEK293NK-SEQ3_RPKM','HEK293-SEQ1_RPKM','HEK293-SEQ2_RPKM','HEK293-SEQ3_RPKM']].values
    rpkm_var = pd.DataFrame(index=['NK.1','NK.2','NK.3','WT.1','WT.2','WT.3']) # annot for variables (cols)

    # convert to AnnData objects (default datatype is 'float32')
    adata_raw = ad.AnnData(raw_X, obs=obs, var=raw_var)
    adata_rpkm = ad.AnnData(rpkm_X, obs=obs, var=rpkm_var)

    return [adata_raw, adata_rpkm]