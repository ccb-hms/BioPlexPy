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