#!/usr/bin/env python

import pandas as pd
import networkx as nx
import numpy as np
import itertools
import random

def bioplex2graph(bp_PPI_df):
    '''
    Convert BioPlex PPIs into a graph.
    
    This function converts representation of BioPlex PPIs into a graph data structure
    representation of BioPlex PPIs in a NetworkX object from NetworkX. 

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

def get_PPI_network_for_complex(bp_PPI_G, Corum_DF, Complex_ID):
    '''
    Retrieve Network of BioPlex (AP-MS) PPIs for a CORUM complex.
    
    This function returns a subgraph of PPIs identified through AP-MS
    between the proteins in a specified CORUM complex.

    Parameters
    ----------
    Network of PPIs : NetworkX graph
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    NetworkX Graph
        A subgraph induced by the proteins in a CORUM complex from the BioPlex network used as input.

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> bp_293t_G = bioplex2graph(bp_293t_df) # (2) Obtain NetworkX graph representation of 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (3) Obtain CORUM complexes
    >>> ING2_bp_293t_G = get_PPI_network_for_complex(bp_293t_G, Corum_DF, 2851) # (4) Get AP-MS interactions as subgraph for a specified protein complex using PPI data
    '''
    # store gene UNIPROT IDs that belong to this complex in a list
    genes_in_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(UniProt IDs)'].values[0].split(';')
    
    # get subgraph induced by the subset of nodes in this CORUM complex
    bp_complex_i_G = bp_PPI_G.subgraph(genes_in_complex_i)
    
    return bp_complex_i_G

def get_prop_edges_in_complex_identfied(bp_PPI_df, Corum_DF, Complex_ID):
    '''
    Calculates proportion of all possible edges identified from BioPlex (AP-MS) PPIs for a CORUM complex.
    
    This function returns the proportion of all possible PPIs identified through AP-MS
    between the proteins in a specified CORUM complex.

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    Float
        The proportion of interactions between all proteins in CORUM complex identified through AP-MS PPI data

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (2) Obtain CORUM complexes
    >>> get_prop_edges_in_complex_identfied(bp_293t_df, Corum_DF, 2851) # (3) Get proportion of interactions identified for a specified CORUM complex using PPI data
    '''
    # store gene symbols that belong to this complex in a list
    genes_in_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(Gene name)'].values[0].split(';')

    # filter BioPlex PPI dataframe to include only interactions where both genes are found in complex
    complex_i_PPI_filter = []
    for symbol_A, symbol_B in zip(bp_PPI_df.SymbolA, bp_PPI_df.SymbolB):

        # check to see if both gene symbols for this interaction are genes in complex
        if (symbol_A in genes_in_complex_i) and (symbol_B in genes_in_complex_i):
            complex_i_PPI_filter.append(True)
        else:
            complex_i_PPI_filter.append(False)

    complex_i_PPI_filter = np.array(complex_i_PPI_filter)
    bp_complex_i_df = bp_PPI_df[complex_i_PPI_filter] # use filter to subset bp PPI dataframe, AP-MS interactions for this CORUM complex
    bp_complex_i_df.reset_index(inplace = True, drop = True) # reset index
    bp_complex_i_df = bp_complex_i_df.loc[:,['SymbolA','SymbolB']] # subset PPI dataframe to the cols we need to construct graph

    # create a graph from the nodes/genes of specified complex
    bp_complex_i_G = nx.Graph()
    bp_complex_i_G.add_nodes_from(genes_in_complex_i)

    # iterate over AP-MS interactions in PPI df and add edges
    for source, target in zip(bp_complex_i_df.SymbolA, bp_complex_i_df.SymbolB):
        bp_complex_i_G.add_edge(source, target)
        
    # create a complete graph from the nodes of complex graph (all possible interactions between proteins)
    bp_complex_i_G_complete = nx.Graph()
    bp_complex_i_G_complete.add_nodes_from(bp_complex_i_G.nodes)
    bp_complex_i_G_complete.add_edges_from(itertools.combinations(bp_complex_i_G.nodes, 2))
    
    # calculate proportion of interactions between proteins in complex identified through AP-MS
    prop_edges_identified = float(len(list(bp_complex_i_G.edges)))/float(len(list(bp_complex_i_G_complete.edges)))
    
    return round(prop_edges_identified, 3) # return proportion of edges ID'd through AP-MS, round to 3 decimal places

def permutation_test_for_CORUM_complex(bp_PPI_df, bp_PPI_G, Corum_DF, Complex_ID):
    '''
    Run permutation test to check for enrichment of BioPlex (AP-MS) PPIs for a given CORUM complex.

    This function returns a p-value after running a permutation test by 1. taking the number of 
    proteins in the specified CORUM complex (N), 2. choosing N random proteins from the Graph generated by 
    all of the PPI data (G), 3. calculating the number of edges in the Subgraph (S) induced by N random proteins
    and storing this value (E_i), 4. repeating steps 1-3 1000 times to create a null distribution, 5. calculating the
    number of edges between N proteins in the CORUM complex (E), 6. returning a p-value by calculating the proportion
    of values [E_1, E_2, ... , E_1000] that are greater than or equal to E.

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame
    Network of PPIs : NetworkX graph
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    Float
        A p-value from a permutation test to check for enrichment of PPIs detected between proteins of CORUM complex

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> bp_293t_G = bioplex2graph(bp_PPI_df) # (2) Obtain NetworkX graph representation of 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (3) Obtain CORUM complexes
    >>> permutation_test_for_CORUM_complex(bp_293t_df, bp_293t_G, Corum_DF, 27) # (4) Calculate p-value to check for enrichment of edges in Arp2/3 protein complex
    '''
    # store gene symbols that belong to this complex in a list
    genes_in_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(Gene name)'].values[0].split(';')

    # filter BioPlex PPI dataframe to include only interactions where both genes are found in complex
    complex_i_PPI_filter = []
    for symbol_A, symbol_B in zip(bp_PPI_df.SymbolA, bp_PPI_df.SymbolB):

        # check to see if both gene symbols for this interaction are genes in complex
        if (symbol_A in genes_in_complex_i) and (symbol_B in genes_in_complex_i):
            complex_i_PPI_filter.append(True)
        else:
            complex_i_PPI_filter.append(False)

    complex_i_PPI_filter = np.array(complex_i_PPI_filter)
    bp_complex_i_df = bp_PPI_df[complex_i_PPI_filter] # use filter to subset bp PPI dataframe, AP-MS interactions for this CORUM complex
    bp_complex_i_df.reset_index(inplace = True, drop = True) # reset index
    bp_complex_i_df = bp_complex_i_df.loc[:,['SymbolA','SymbolB']] # subset PPI dataframe to the cols we need to construct graph

    # create a graph from the nodes/genes of specified complex
    bp_complex_i_G = nx.Graph()
    bp_complex_i_G.add_nodes_from(genes_in_complex_i)

    # iterate over AP-MS interactions in PPI df and add edges
    for source, target in zip(bp_complex_i_df.SymbolA, bp_complex_i_df.SymbolB):
        bp_complex_i_G.add_edge(source, target)

    # number of edges detected between proteins in this CORUM complex among PPI data
    num_edges_identified_CORUM_complex = float(len(list(bp_complex_i_G.edges)))

    # number of genes in CORUM complex (N genes)
    num_genes_in_CORUM_complex = len(list(bp_complex_i_G.nodes))

    # list of nodes in large network generated from PPI data
    nodes_in_overall_PPI_network = list(bp_PPI_G.nodes)

    # list that will store number of edges detected in each subgraph
    num_edges_random_subgraphs = []

    # iterate through 1000 random subgraphs induced by N nodes
    for S_i in np.arange(0, 1000):

        # choose N genes at random without replacement
        N_rando_nodes_from_PPI_network = random.sample(nodes_in_overall_PPI_network, num_genes_in_CORUM_complex)

        # get subgraph induced by random subset of nodes
        bp_PPI_S = bp_PPI_G.subgraph(N_rando_nodes_from_PPI_network)

        # calculate the number of edges detected within subgraph induced by random nodes
        num_edges_S = float(len(list(bp_PPI_S.edges)))

        # store in list that contains permutations
        num_edges_random_subgraphs.append(num_edges_S)

    # convert list to numpy array
    num_edges_random_subgraphs = np.array(num_edges_random_subgraphs)

    # calculate proportion of subgraphs that had more edges than edges detected in CORUM complex (p-val from permutation test)
    p_val = float(np.sum(num_edges_random_subgraphs >= num_edges_identified_CORUM_complex)) / 1000.0

    return p_val