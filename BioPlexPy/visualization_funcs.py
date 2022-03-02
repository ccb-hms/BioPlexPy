#!/usr/bin/env python

import itertools
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def display_PPI_network_for_complex(ax, bp_PPI_df, Corum_DF, Complex_ID, node_size, edge_width, node_font_size = 10, bait_node_color='xkcd:red', prey_node_color='xkcd:rose pink', AP_MS_edge_color='xkcd:red', node_pos=False):
    '''
    Display network of BioPlex PPIs for a CORUM complex.
    
    This function displays a complete network in which nodes represent 
    the proteins in a specified CORUM complex and edges represent 
    BioPlex PPIs using NetworkX. Edges detected through AP-MS are colored darker.

    Parameters
    ----------
    ax object to draw on: Matplotlib Axes
    DataFrame of PPIs : Pandas DataFrame
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int
    Size of Nodes in Network: int
    Width of Edges in Network: float
    Size of font for Node Labels: int (optional)
    Color of Nodes targeted as baits: str (optional)
    Color of Nodes detected as preys only: str (optional)
    Color of Edges observed via AP-MS from PPI data: str (optional)
    Networkx Position of Nodes: dict (optional)

    Returns
    -------
    Node Positions
        Dictionary of Node Positions in NetworkX layout

    Examples
    --------
    >>> bp_PPI_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (2) Obtain CORUM complexes
    >>> ING2_node_layout = display_PPI_network_for_complex(bp_PPI_df, Corum_DF, 2851, 2300, 3.5) # (3) Visualize network for specified protein complex using PPI data (ING2 complex ID: 2851)
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
    bp_complex_i_df = bp_PPI_df[complex_i_PPI_filter] # use filter to subset bp PPI dataframe
    bp_complex_i_df.reset_index(inplace = True, drop = True) # reset index
    bp_complex_i_df = bp_complex_i_df.loc[:,['SymbolA','SymbolB']] # subset PPI dataframe to the cols we need to construct graph

    # create a graph from the nodes/genes of specified complex
    bp_complex_i_G = nx.Graph()
    bp_complex_i_G.add_nodes_from(genes_in_complex_i)

    # iterate over AP-MS interactions in PPI df and add edges
    for source, target in zip(bp_complex_i_df.SymbolA, bp_complex_i_df.SymbolB):
        bp_complex_i_G.add_edge(source, target)

    # get a list of genes that were identifed as "baits" and "preys" for coloring nodes
    bp_complex_i_baits = list(set(bp_complex_i_df.SymbolA))
    bp_complex_i_preys = list(set(bp_complex_i_df.SymbolB))

    #labels = nx.get_node_attributes(bp_complex_i_G, 'symbol') # DEFAULT: labels will be the gene symbols, pull from attribute for each node
    labels = dict([(key,val) for key, val in zip(list(bp_complex_i_G.nodes), list(bp_complex_i_G.nodes))])

    # color nodes according to whether they were present among "baits" & "preys", just "baits", just "preys" or not detected in PPI data for this complex
    node_color_map = []
    for node_i_uniprot in bp_complex_i_G.nodes:

        node_i_gene_symbol = labels[node_i_uniprot] # get gene symbol

        # gene is present among baits & preys for the PPIs detected in this complex
        if (node_i_gene_symbol in bp_complex_i_baits) and (node_i_gene_symbol in bp_complex_i_preys):
            node_color_map.append(bait_node_color)

        # gene is present among baits but NOT preys for the PPIs detected in this complex
        elif (node_i_gene_symbol in bp_complex_i_baits) and (node_i_gene_symbol not in bp_complex_i_preys):
            node_color_map.append(bait_node_color)

        # gene is NOT present among baits but is present among preys for the PPIs detected in this complex
        elif (node_i_gene_symbol not in bp_complex_i_baits) and (node_i_gene_symbol in bp_complex_i_preys):
            node_color_map.append(prey_node_color)

        # gene is NOT present among baits and is NOT present among preys for the PPIs detected in this complex
        elif (node_i_gene_symbol not in bp_complex_i_baits) and (node_i_gene_symbol not in bp_complex_i_preys):
            node_color_map.append('0.7')

    # create a complete graph from the nodes of complex graph to add in all "background" edges (edges detected with AP-MS will be colored over)
    # position will be the same for both graphs since nodes are the same
    bp_complex_i_G_complete = nx.Graph()
    bp_complex_i_G_complete.add_nodes_from(bp_complex_i_G.nodes)
    bp_complex_i_G_complete.add_edges_from(itertools.combinations(bp_complex_i_G.nodes, 2))

    ############################################################################
    #### optional argument "node_pos" used here
    ############################################################################
    # check to see if node position object has been fed as an argument
    if node_pos == False:
        pos = nx.circular_layout(bp_complex_i_G)  # DEFAULT: set position of nodes w/ circular layout
    else:
        pos = node_pos # use node positions fed into function
    ############################################################################

    # construct edges for COMPLETE graph for "background" edges
    edges_complete = nx.draw_networkx_edges(bp_complex_i_G_complete, pos, width = edge_width, alpha = 0.25, ax = ax)
    edges_complete.set_edgecolor("xkcd:grey")

    # construct edges
    edges = nx.draw_networkx_edges(bp_complex_i_G, pos, width = edge_width, ax = ax)
    edges.set_edgecolor(AP_MS_edge_color)

    # construct nodes
    nodes = nx.draw_networkx_nodes(bp_complex_i_G, pos, node_size = node_size, node_color = node_color_map, ax = ax)
    nodes.set_edgecolor("xkcd:black")
    nodes.set_linewidth(1.5)
    nx.draw_networkx_labels(bp_complex_i_G, pos, labels = labels, font_size = node_font_size, font_weight = 'bold', font_color = 'xkcd:white', ax = ax)

    # return node position layout
    return pos