#!/usr/bin/env python

import numpy as np
import itertools
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
    >>> fig, ax = plt.subplots() # (3) create figure and axis objects to draw on
    >>> ING2_node_layout = display_PPI_network_for_complex(ax, bp_PPI_df, Corum_DF, 2851, 2300, 3.5) # (4) Visualize network for specified protein complex using PPI data (ING2 complex ID: 2851)
    '''
    # store uniprot IDs & gene symbols that belong to this complex in a list
    genes_in_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(UniProt IDs)'].values[0].split(';') # Uniprot
    gene_symbols_in_complex_i = Corum_DF[Corum_DF.ComplexID == Complex_ID].loc[:,'subunits(Gene name)'].values[0].split(';') # Symbol

    # filter BioPlex PPI dataframe to include only interactions where both genes are found in complex
    complex_i_PPI_filter = []
    for uniprot_A, uniprot_B in zip(bp_PPI_df.UniprotA, bp_PPI_df.UniprotB):
        
        # check for isoform IDs and adjust
        if '-' in uniprot_A:
            uniprot_A = uniprot_A.split('-')[0]
        if '-' in uniprot_B:
            uniprot_B = uniprot_B.split('-')[0]

        # check to see if both gene symbols for this interaction are genes in complex
        if (uniprot_A in genes_in_complex_i) and (uniprot_B in genes_in_complex_i):
            complex_i_PPI_filter.append(True)
        else:
            complex_i_PPI_filter.append(False)

    complex_i_PPI_filter = np.array(complex_i_PPI_filter)
    bp_complex_i_df = bp_PPI_df[complex_i_PPI_filter] # use filter to subset bp PPI dataframe
    bp_complex_i_df.reset_index(inplace = True, drop = True) # reset index

    # reconstruct UniprotA/UniprotB columns without '-' isoform id
    UniprotA_new = []
    UniprotB_new = []
    for UniprotA, UniprotB in zip(bp_complex_i_df.UniprotA, bp_complex_i_df.UniprotB):

        if '-' in UniprotA:
            UniprotA_new.append(UniprotA.split('-')[0])
        else:
            UniprotA_new.append(UniprotA)

        if '-' in UniprotB:
            UniprotB_new.append(UniprotB.split('-')[0])
        else:
            UniprotB_new.append(UniprotB)
            
    # update columns for Uniprot source & Uniprot target to exclude isoform '-' ID
    bp_complex_i_df.loc[:,'UniprotA'] = UniprotA_new
    bp_complex_i_df.loc[:,'UniprotB'] = UniprotB_new
    
    # subset PPI dataframe to the cols we need to construct graph
    bp_complex_i_df = bp_complex_i_df.loc[:,['UniprotA','UniprotB','SymbolA','SymbolB']]
    
    # create a graph from the nodes/genes of specified complex
    bp_complex_i_G = nx.Graph()
    bp_complex_i_G.add_nodes_from(genes_in_complex_i)
    
    # iterate over AP-MS interactions in PPI df and add edges
    for source, target in zip(bp_complex_i_df.UniprotA, bp_complex_i_df.UniprotB):
        bp_complex_i_G.add_edge(source, target)
        
    # get mapping uniprot -> symbol & store as node attribute from CORUM complex data
    uniprot_symbol_dict = dict([key,val] for key, val in zip(genes_in_complex_i, gene_symbols_in_complex_i))
    for node_i in bp_complex_i_G.nodes():
        bp_complex_i_G.nodes[node_i]["symbol"] = uniprot_symbol_dict[node_i]

    # get a list of genes that were identifed as "baits" and "preys" for coloring nodes
    bp_complex_i_baits = list(set(bp_complex_i_df.UniprotA))
    bp_complex_i_preys = list(set(bp_complex_i_df.UniprotB))

    labels = dict([(key,val) for key, val in zip(list(bp_complex_i_G.nodes), [bp_complex_i_G.nodes[node_i]['symbol'] for node_i in bp_complex_i_G.nodes])])

    # color nodes according to whether they were present among "baits" & "preys", just "baits", just "preys" or not detected in PPI data for this complex
    node_color_map = []
    for node_i_uniprot in bp_complex_i_G.nodes:

        # gene is present among baits & preys for the PPIs detected in this complex
        if (node_i_uniprot in bp_complex_i_baits) and (node_i_uniprot in bp_complex_i_preys):
            node_color_map.append(bait_node_color)

        # gene is present among baits but NOT preys for the PPIs detected in this complex
        elif (node_i_uniprot in bp_complex_i_baits) and (node_i_uniprot not in bp_complex_i_preys):
            node_color_map.append(bait_node_color)

        # gene is NOT present among baits but is present among preys for the PPIs detected in this complex
        elif (node_i_uniprot not in bp_complex_i_baits) and (node_i_uniprot in bp_complex_i_preys):
            node_color_map.append(prey_node_color)

        # gene is NOT present among baits and is NOT present among preys for the PPIs detected in this complex
        elif (node_i_uniprot not in bp_complex_i_baits) and (node_i_uniprot not in bp_complex_i_preys):
            node_color_map.append('0.7')

    # create a complete graph from the nodes of complex graph to add in all "background" edges (edges detected with AP-MS will be colored over)
    # position will be the same for both graphs since nodes are the same
    bp_complex_i_G_complete = nx.Graph()
    bp_complex_i_G_complete.add_nodes_from(bp_complex_i_G.nodes)
    bp_complex_i_G_complete.add_edges_from(itertools.combinations(bp_complex_i_G.nodes, 2))

    # optional argument "node_pos" used here
    # check to see if node position object has been fed as an argument
    if node_pos == False:
        pos = nx.circular_layout(bp_complex_i_G)  # DEFAULT: set position of nodes w/ circular layout
    else:
        pos = node_pos # use node positions fed into function

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

def display_PDB_network_for_complex(ax, interacting_UniProt_IDs, node_size, edge_width, node_font_size=10):
    '''
    Display network of interacting chains.
    
    This function displays a complete network in which nodes represent the 
    proteins in a specified PDB structure, and edges represent chains in that
    structure, using NetworkX. Edges that are classified as interacting
    (are < 6 angstroms apart) are colored black.

    Parameters
    ----------
    ax object to draw on: Matplotlib Axes
    List of Interacting Chains: list
    Size of Nodes in Network: int
    Width of Edges in Network: float
    Size of font for Node Labels: int (optional)

    Returns
    -------
    Node Positions
        Dictionary of Node Positions in NetworkX layout
    Interacting Network Edges
        List of Edges for Interacting Nodes
    Number of Network Edges
        Float of the Number of Possible Interacting Edges

    Examples
    --------
    >>> interacting_chains_list = get_interacting_chains_from_PDB('6YW7', '/n/data1/hms/ccb/lab/projects/bioplex/BioPlexPy/protein_function_testing') # (1) Obtain list of interacting chains from 6YW7 structure
    >>> chain_to_UniProt_mapping_dict = list_uniprot_pdb_mappings('6YW7') # (2) Obtain a mapping of PDB ID 6YW7 chains to UniProt IDs
    >>> interacting_UniProt_IDs = PDB_chains_to_uniprot(interacting_chains_list, chain_to_UniProt_mapping_dict) # (3) Obtain list of interacting chains from 6YW7 structure using UniProt IDs
    >>> fig, ax = plt.subplots() # (4) create figure and axis objects to draw on
    >>> node_layout_pdb, edges_list_pdb, num_possible_edges_pdb = display_PDB_network_for_complex(ax, interacting_UniProt_IDs, 2300, 3.5) # (5) Visualize interacting chains using Uniprot IDs
    '''
    # create connected graph from all uniprot IDs
    chain_uniprot_IDs = [chain_uniprot_i[0] for chain_uniprot_i in chain_to_UniProt_mapping_dict.values()]

    # create a graph from the nodes/genes of this complex structure
    pdb_structure_i_G = nx.Graph()
    pdb_structure_i_G.add_nodes_from(chain_uniprot_IDs)

    # iterate over physical interactions and add edges
    for chain_pair_interacting in interacting_UniProt_IDs:
        pdb_structure_i_G.add_edge(chain_pair_interacting[0], chain_pair_interacting[1])

    # create a complete graph from the nodes of complex graph to add in all "background" edges (edges detected with AP-MS will be colored over)
    # position will be the same for both graphs since nodes are the same
    pdb_structure_i_G_complete = nx.Graph()
    pdb_structure_i_G_complete.add_nodes_from(pdb_structure_i_G.nodes)
    pdb_structure_i_G_complete.add_edges_from(itertools.combinations(pdb_structure_i_G.nodes, 2))

    # set position of nodes w/ circular layout
    pos = nx.circular_layout(pdb_structure_i_G)

    # construct edges for COMPLETE graph for "background" edges
    edges_complete = nx.draw_networkx_edges(pdb_structure_i_G_complete, pos, width = edge_width, alpha = 0.25, ax = ax)
    edges_complete.set_edgecolor("xkcd:grey")

    # construct edges
    edges = nx.draw_networkx_edges(pdb_structure_i_G, pos, width = edge_width, ax = ax)
    edges.set_edgecolor('xkcd:blue')

    # construct nodes
    nodes = nx.draw_networkx_nodes(pdb_structure_i_G, pos, node_size = node_size, node_color = 'xkcd:black', ax = ax)
    nodes.set_edgecolor("xkcd:black")
    nodes.set_linewidth(1.5)
    nx.draw_networkx_labels(pdb_structure_i_G, pos, font_size = node_font_size, font_weight = 'bold', font_color = 'xkcd:white', ax = ax)

    # return node position layout, list of edges detected, number of possible edges
    return [pos, pdb_structure_i_G.edges, float(len(pdb_structure_i_G_complete.edges))]

def display_PPI_network_match_PDB(ax, interacting_UniProt_IDs, bp_PPI_df, node_pos, node_size, edge_width, node_font_size=10, bait_node_color='xkcd:red', prey_node_color='xkcd:rose pink', AP_MS_edge_color='xkcd:red'):
    '''
    Display network of BioPlex PPIs for a set of interacting UniProt IDs.
    
    This function displays a complete network in which nodes represent the 
    proteins in a specified PDB structure, and edges represent chains in that
    structure, using NetworkX. Edges that are classified as interacting from
    BioPlex PPI data (detected through AP-MS) are colored darker.

    Parameters
    ----------
    ax object to draw on: Matplotlib Axes
    List of Interacting Chains: list
    DataFrame of PPIs : Pandas DataFrame
    Networkx Position of Nodes: dict
    Size of Nodes in Network: int
    Width of Edges in Network: float
    Size of font for Node Labels: int (optional)
    Color of Nodes targeted as baits: str (optional)
    Color of Nodes detected as preys only: str (optional)
    Color of Edges observed via AP-MS from PPI data: str (optional)

    Returns
    -------
    Interacting Network Edges
        List of Edges for Interacting Nodes
    Number of Network Edges
        Float of the Number of Possible Interacting Edges

    Examples
    --------
    >>> interacting_chains_list = get_interacting_chains_from_PDB('6YW7', '/n/data1/hms/ccb/lab/projects/bioplex/BioPlexPy/protein_function_testing') # (1) Obtain list of interacting chains from 6YW7 structure
    >>> chain_to_UniProt_mapping_dict = list_uniprot_pdb_mappings('6YW7') # (2) Obtain a mapping of PDB ID 6YW7 chains to UniProt IDs
    >>> interacting_UniProt_IDs = PDB_chains_to_uniprot(interacting_chains_list, chain_to_UniProt_mapping_dict) # (3) Obtain list of interacting chains from 6YW7 structure using UniProt IDs
    >>> fig, ax1 = plt.subplots() # (4) create figure and axis objects to draw on
    >>> node_layout_pdb, edges_list_pdb, num_possible_edges_pdb = display_PDB_network_for_complex(ax1, interacting_UniProt_IDs, 2300, 3.5) # (5) Visualize interacting chains using Uniprot IDs
    >>> bp_PPI_df = getBioPlex('293T', '3.0') # (6) Get BioPlex PPI data
    >>> fig, ax2 = plt.subplots() # (7) create figure and axis objects to draw on
    >>> edges_list_bp, num_possible_edges_bp = display_PPI_network_match_PDB(ax2, interacting_UniProt_IDs, bp_PPI_df, node_layout_pdb, 2300, 3.5) # (8) Visualize BioPlex PPI interactions using layout from interacting chains
    '''
    # create connected graph from all uniprot IDs
    chain_uniprot_IDs = [chain_uniprot_i[0] for chain_uniprot_i in chain_to_UniProt_mapping_dict.values()]
    
    # filter BioPlex PPI dataframe to include only interactions where both genes are found in complex
    structure_uniprots_i_PPI_filter = []
    for uniprot_A, uniprot_B in zip(bp_PPI_df.UniprotA, bp_PPI_df.UniprotB):

        # check for isoform IDs and adjust
        if '-' in uniprot_A:
            uniprot_A = uniprot_A.split('-')[0]
        if '-' in uniprot_B:
            uniprot_B = uniprot_B.split('-')[0]

        # check to see if both gene symbols for this interaction are genes in complex
        if (uniprot_A in chain_uniprot_IDs) and (uniprot_B in chain_uniprot_IDs):
            structure_uniprots_i_PPI_filter.append(True)
        else:
            structure_uniprots_i_PPI_filter.append(False)

    structure_uniprots_i_PPI_filter = np.array(structure_uniprots_i_PPI_filter)
    bp_structure_i_df = bp_PPI_df[structure_uniprots_i_PPI_filter] # use filter to subset bp PPI dataframe
    bp_structure_i_df.reset_index(inplace = True, drop = True) # reset index

    # reconstruct UniprotA/UniprotB columns without '-' isoform id
    UniprotA_new = []
    UniprotB_new = []
    for UniprotA, UniprotB in zip(bp_structure_i_df.UniprotA, bp_structure_i_df.UniprotB):

        if '-' in UniprotA:
            UniprotA_new.append(UniprotA.split('-')[0])
        else:
            UniprotA_new.append(UniprotA)

        if '-' in UniprotB:
            UniprotB_new.append(UniprotB.split('-')[0])
        else:
            UniprotB_new.append(UniprotB)

    # update columns for Uniprot source & Uniprot target to exclude isoform '-' ID
    bp_structure_i_df.loc[:,'UniprotA'] = UniprotA_new
    bp_structure_i_df.loc[:,'UniprotB'] = UniprotB_new

    # subset PPI dataframe to the cols we need to construct graph
    bp_structure_i_df = bp_structure_i_df.loc[:,['UniprotA','UniprotB','SymbolA','SymbolB']]

    # create a graph from the nodes/genes of specified complex
    bp_structure_i_G = nx.Graph()
    bp_structure_i_G.add_nodes_from(chain_uniprot_IDs)

    # iterate over AP-MS interactions in PPI df and add edges
    for source, target in zip(bp_structure_i_df.UniprotA, bp_structure_i_df.UniprotB):
        bp_structure_i_G.add_edge(source, target)

    # get a list of genes that were identifed as "baits" and "preys" for coloring nodes
    bp_structure_i_baits = list(set(bp_structure_i_df.UniprotA))
    bp_structure_i_preys = list(set(bp_structure_i_df.UniprotB))

    # color nodes according to whether they were present among "baits" & "preys", just "baits", just "preys" or not detected in PPI data for this structure
    node_color_map = []
    for node_i_uniprot in bp_structure_i_G.nodes:

        # gene is present among baits & preys for the PPIs detected in this structure
        if (node_i_uniprot in bp_structure_i_baits) and (node_i_uniprot in bp_structure_i_preys):
            node_color_map.append(bait_node_color)

        # gene is present among baits but NOT preys for the PPIs detected in this structure
        elif (node_i_uniprot in bp_structure_i_baits) and (node_i_uniprot not in bp_structure_i_preys):
            node_color_map.append(bait_node_color)

        # gene is NOT present among baits but is present among preys for the PPIs detected in this structure
        elif (node_i_uniprot not in bp_structure_i_baits) and (node_i_uniprot in bp_structure_i_preys):
            node_color_map.append(prey_node_color)

        # gene is NOT present among baits and is NOT present among preys for the PPIs detected in this complex
        elif (node_i_uniprot not in bp_structure_i_baits) and (node_i_uniprot not in bp_structure_i_preys):
            node_color_map.append('0.7')
            
    # can use the complete graph from the pdb direct interaction graph since they have the same nodes
    # create a complete graph from the nodes of complex graph to add in all "background" edges (edges detected with AP-MS will be colored over)
    # position will be the same for both graphs since nodes are the same
    bp_structure_i_G_complete = nx.Graph()
    bp_structure_i_G_complete.add_nodes_from(bp_structure_i_G.nodes)
    bp_structure_i_G_complete.add_edges_from(itertools.combinations(bp_structure_i_G.nodes, 2))

    # construct edges for COMPLETE graph for "background" edges
    edges_complete = nx.draw_networkx_edges(bp_structure_i_G_complete, node_pos, width = edge_width, alpha = 0.25, ax = ax)
    edges_complete.set_edgecolor("xkcd:grey")

    # construct edges
    edges = nx.draw_networkx_edges(bp_structure_i_G, node_pos, width = edge_width, ax = ax)
    edges.set_edgecolor(AP_MS_edge_color)

    # construct nodes
    nodes = nx.draw_networkx_nodes(bp_structure_i_G, node_pos, node_size = node_size, node_color = node_color_map, ax = ax)
    nodes.set_edgecolor("xkcd:black")
    nodes.set_linewidth(1.5)
    nx.draw_networkx_labels(bp_structure_i_G, node_pos, font_size = node_font_size, font_weight = 'bold', font_color = 'xkcd:white', ax = ax)

    # return node position layout, list of edges detected, number of possible edges
    return [bp_structure_i_G.edges, float(len(bp_structure_i_G_complete.edges))]