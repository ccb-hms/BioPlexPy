#!/usr/bin/env python

import pandas as pd
import networkx as nx

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
