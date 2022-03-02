#!/usr/bin/env python

import pandas as pd
import networkx as nx
import numpy as np
import itertools

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

def get_PPI_network_edges_for_complex(bp_PPI_df, Corum_DF, Complex_ID):
    '''
    Retrieve DataFrame of BioPlex (AP-MS) PPIs for a CORUM complex.
    
    This function returns a dataframe of PPIs identified through AP-MS
    between the proteins in a specified CORUM complex.

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    Pandas DataFrame
        A dataframe with each row corresponding to a PPI interaction detected via AP-MS for specified complex

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (2) Obtain CORUM complexes
    >>> get_PPI_network_edges_for_complex(bp_293t_df, Corum_DF, 2851) # (3) Get AP-MS interactions for a specified protein complex using PPI data
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
    bp_complex_i_df = bp_complex_i_df.loc[:,['SymbolA', 'SymbolB', 'pW', 'pNI', 'pInt']] # drop columns not needed
    
    return bp_complex_i_df

def get_prop_edges_in_complex_identfied(bp_PPI_df, Corum_DF, Complex_ID):
    '''
    Retrieve DataFrame of BioPlex (AP-MS) PPIs for a CORUM complex.
    
    This function returns a dataframe of PPIs identified through AP-MS
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