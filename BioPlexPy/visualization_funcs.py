#!/usr/bin/env python

import itertools
import numpy as np
import matplotlib.pyplot as plt

def display_PPI_network_for_complex(bp_PPI_df, Corum_DF, ComplexName_i, node_size, edge_width, bait_node_color, prey_node_color, AP_MS_edge_color, fig_out_path, node_pos = False, node_labels = False):
    '''
    Display network of BioPlex PPIs for a CORUM complex.
    
    This function displays a complete network in which nodes represent 
    the proteins in a specified CORUM complex and edges represent 
    BioPlex PPIs using NetworkX. Edges detected through AP-MS are colored darker.

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame
    DataFrame of CORUM complexes : Pandas DataFrame
    Name of Corum Complex: str
    Size of Nodes in Network: int
    Width of Edges in Network: float
    Color of Nodes targeted as baits: str
    Color of Nodes detected as preys only: str
    Color of Edges observed via AP-MS from PPI data: str
    Path to save figure: str
    Networkx Position of Nodes: dict (optional)
    Node Labels: dict (optional)

    Returns
    -------
    PNG file
        Displays and outputs a PNG file of network.
    Node Positions
        Dictionary of Node Positions in NetworkX layout
    Node Labels
        Dictionary of Labels for Nodes used in plotting

    Examples
    --------
    >>> bp_293t_df = getBioPlex('293T', '3.0') # (1) Obtain the latest version of the 293T PPI network
    >>> Corum_DF = getCorum('core', 'Human') # (2) Obtain CORUM complexes
    >>> ING2_node_layout, ING2_node_labels = display_PPI_network_for_complex(bp_PPI_df, Corum_DF, 'ING2 complex', 2300, 3.5, /n/shared_db/ccb/bioplex/BioPlexPy_testing/figures/network_293T_3.0_ING2-complex.png') # (3) Visualize network for specified protein complex using PPI data
    '''
    # store gene symbols that belong to this complex in a list
    genes_in_complex_i = Corum_DF[Corum_DF.ComplexName == ComplexName_i].loc[:,'subunits(Gene name)'].values[0].split(';')

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
    bp_complex_i_G = bioplex2graph(bp_complex_i_df) # turn subsetted bp PPI df into a graph
    bp_complex_i_G = bp_complex_i_G.to_undirected() # convert to undirected graph for displaying

    # get a list of genes that were identifed as "baits" and "preys" for coloring nodes
    bp_complex_i_baits = list(set(bp_complex_i_df.SymbolA))
    bp_complex_i_preys = list(set(bp_complex_i_df.SymbolB))

    ############################################################################
    #### optional argument "node_labels" used here
    ############################################################################
    # check to see if node labels have been fed as an argument
    if node_labels == False:
        labels = nx.get_node_attributes(bp_complex_i_G, 'symbol') # DEFAULT: labels will be the gene symbols, pull from attribute for each node

        # add genes in CORUM complex that were not detected as baits or preys as "dummy" nodes (for visualization purposes only)
        dummy_node = 1
        for gene_i in genes_in_complex_i:
            if gene_i not in list(labels.values()): # gene already not a node
                bp_complex_i_G.add_nodes_from([(f'dummy_{dummy_node}', {"symbol": gene_i})])
                labels[f'dummy_{dummy_node}'] = gene_i
                dummy_node += 1
                
    # use labels passed as argument to keep things consistent between network visualizations
    else:
        labels = node_labels # use node labels fed into function
        labels_r = {value : key for (key, value) in labels.items()} # reverse dictionary to add nodes to graph
        
        # add genes in CORUM complex that were not detected as baits or preys as "dummy" nodes (for visualization purposes only)
        for gene_i in genes_in_complex_i:
            if gene_i not in list(bp_complex_i_G.nodes()): # gene not already a node
                bp_complex_i_G.add_nodes_from([(labels_r[gene_i], {"symbol": gene_i})])
    ############################################################################

    # color nodes according to whether they were present among "baits" & "preys", just "baits", or just "preys" for this complex
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
            
        # gene is NOT present among baits and is NOT present among preys for the PPIs detected in this complex ("dummy" nodes)
        elif (node_i_gene_symbol not in bp_complex_i_baits) and (node_i_gene_symbol not in bp_complex_i_preys):
            node_color_map.append('0.7')

    # create a complete graph from the nodes of complex graph to add in all "background" edges (edges detected with AP-MS will be colored over)
    # position will be the same for both graphs since nodes are the same
    bp_complex_i_G_complete = nx.Graph()
    bp_complex_i_G_complete.add_nodes_from(bp_complex_i_G.nodes)
    bp_complex_i_G_complete.add_edges_from(itertools.combinations(bp_complex_i_G.nodes, 2))

    # create figure instance
    fig, ax = plt.subplots()
    
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
    edges_complete = nx.draw_networkx_edges(bp_complex_i_G_complete, pos, width = edge_width, alpha = 0.25)
    edges_complete.set_edgecolor("xkcd:grey")

    # construct edges
    edges = nx.draw_networkx_edges(bp_complex_i_G, pos, width = edge_width)
    edges.set_edgecolor(AP_MS_edge_color)

    # construct nodes
    nodes = nx.draw_networkx_nodes(bp_complex_i_G, pos, node_size = node_size, node_color = node_color_map)
    nodes.set_edgecolor("xkcd:black")
    nodes.set_linewidth(1.5)
    nx.draw_networkx_labels(bp_complex_i_G, pos, labels = labels, font_size = 10, font_weight = 'bold', font_color = 'xkcd:white')

    fig = plt.gcf()
    fig.set_size_inches(7.5, 7.5)
    fig.tight_layout()
    
    # save figure as PNG
    plt.savefig(fig_out_path, bbox_inches='tight', dpi = 300 , transparent = True)
    plt.show()
    
    # return node position layout & labels for nodes
    return [pos, labels]