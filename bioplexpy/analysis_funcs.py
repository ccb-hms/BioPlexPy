#!/usr/bin/env python

import itertools
import random
import re
from collections import Counter

import networkx as nx
import numpy as np
import pandas as pd
import requests
from Bio.PDB import *
from scipy.spatial.distance import cdist


def bioplex2graph(bp_PPI_df):
    '''
    Convert BioPlex PPIs into a graph.
    
    This function converts representation of BioPlex PPIs into a 
    graph data structure representation of BioPlex PPIs in a NetworkX 
    object from NetworkX. 

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame

    Returns
    -------
    NetworkX graph
        A NetworkX graph with Nodes = Uniprot Gene Symbols 
        and Edges = interactions.

    Examples
    --------
    # (1) Obtain the latest version of the 293T PPI network
    # (2) Turn the data into a graph
    >>> bp_293t_df = getBioPlex('293T', '3.0')
    >>> bp_293t_G = bioplex2graph(bp_293t_df) 
    >>> type(bp_293t_G)
    <class 'networkx.classes.digraph.DiGraph'>
    >>> len(bp_293t_G.edges)
    115868
    >>> len(bp_293t_G.nodes)
    13689
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

    # update columns for Uniprot source 
    # & Uniprot target to exclude isoform '-' ID
    bp_PPI_df.loc[:,'UniprotA'] = UniprotA_new
    bp_PPI_df.loc[:,'UniprotB'] = UniprotB_new
    
    # construct graph from BioPlex PPI data
    bp_G = nx.DiGraph()
    for source, target, pW, pNI, pInt in zip(bp_PPI_df.UniprotA, 
            bp_PPI_df.UniprotB, bp_PPI_df.pW, bp_PPI_df.pNI, bp_PPI_df.pInt):
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
        
    # get set of baits & store a bait boolean as node attribute 
    # if node is a bait True & False otherwise
    bp_i_baits = set(bp_PPI_df.UniprotA)
    for node_i in bp_G.nodes():
        if node_i in bp_i_baits:
            bp_G.nodes[node_i]["bait"] = True
        else:
            bp_G.nodes[node_i]["bait"] = False
    
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
        A subgraph induced by the proteins in a CORUM complex 
        from the BioPlex network used as input.

    Examples
    --------
    # (1) Obtain the latest version of the 293T PPI network
    # (2) Obtain NetworkX graph representation of 293T PPI network
    # (3) Obtain CORUM complexes
    # (4) Get AP-MS interactions as subgraph for a specified protein complex using PPI data
    >>> bp_293t_df = getBioPlex('293T', '3.0')
    >>> bp_293t_G = bioplex2graph(bp_293t_df)
    >>> Corum_DF = getCorum()
    >>> ING2_bp_293t_G = get_PPI_network_for_complex(bp_293t_G, Corum_DF, 2851)
    >>> type(ING2_bp_293t_G)
    <class 'networkx.classes.digraph.DiGraph'>
    >>> len(ING2_bp_293t_G)
    12
    '''
    # store gene UNIPROT IDs that belong to this complex in a list
    genes_in_complex_i = (Corum_DF[Corum_DF.complex_id == Complex_ID].loc[:,
                                'subunits_uniprot_id'].values[0].split(';'))
    
    # get subgraph induced by the subset of nodes in this CORUM complex
    bp_complex_i_G = bp_PPI_G.subgraph(genes_in_complex_i)
    
    return bp_complex_i_G

def get_DataFrame_from_PPI_network(bp_PPI_G):
    '''
    Convert Network of BioPlex (AP-MS) PPIs into DataFrame of BioPlex 
    interaction Network.
    
    This function returns a DataFrame of PPIs (identified through AP-MS) 
    represented as a graph.

    Parameters
    ----------
    Network of PPIs : NetworkX graph

    Returns
    -------
    Pandas DataFrame
        A DataFrame of edges (AP-MS interactions) from a network.

    Examples
    --------
    # (1) Obtain the latest version of the 293T PPI network
    # (2) Obtain NetworkX graph representation of 293T PPI network
    # (3) Obtain CORUM complexes
    # (4) Get AP-MS interactions as subgraph for a specified protein complex using PPI data row corresponding to an edge
    # (5) Convert ING2 AP-MS network into DataFrame w/ each 
    # order of rows in dataframe?
    >>> bp_293t_df = getBioPlex('293T', '3.0')
    >>> bp_293t_G = bioplex2graph(bp_293t_df)
    >>> Corum_DF = getCorum()
    >>> ING2_bp_293t_G = get_PPI_network_for_complex(bp_293t_G, Corum_DF, 2851)
    >>> ING2_bp_293t_df = get_DataFrame_from_PPI_network(ING2_bp_293t_G)
    >>> type(ING2_bp_293t_df)
    <class 'pandas.core.frame.DataFrame'>
    >>> ING2_bp_293t_df.shape
    (51, 7)
    >>> ING2_bp_293t_df.columns
    Index(['UniprotA', 'UniprotB', 'SymbolA', 'SymbolB', 'pW', 'pNI', 'pInt'], dtype='object')
    '''
    # get list of edges in network
    PPI_edge_list = list(bp_PPI_G.edges)

    # make node_A and node_B columns for both UNIPROT & SYMBOLS
    uniprotA_list = []
    uniprotB_list = []
    symbolA_list = []
    symbolB_list = []

    # make columns for calcs detected for each edge
    pW_list = []
    pNI_list = []
    pInt_list = []

    # iterate through each edge to store data for each row of DataFrame
    for edge_i in PPI_edge_list:

        nodeA, nodeB = edge_i

        # nodes are labeled with UNIPROT
        uniprotA_list.append(nodeA)
        uniprotB_list.append(nodeB)

        # get gene SYMBOL
        symbolA_list.append(bp_PPI_G.nodes[nodeA]['symbol'])
        symbolB_list.append(bp_PPI_G.nodes[nodeB]['symbol'])

        # get AP-MS calculations
        pW_list.append(bp_PPI_G.edges[edge_i]['pW'])
        pNI_list.append(bp_PPI_G.edges[edge_i]['pNI'])
        pInt_list.append(bp_PPI_G.edges[edge_i]['pInt'])

    # convert lists into cols of DataFrame
    bp_complex_i_df = pd.DataFrame()
    bp_complex_i_df.loc[:,'UniprotA'] = uniprotA_list
    bp_complex_i_df.loc[:,'UniprotB'] = uniprotB_list
    bp_complex_i_df.loc[:,'SymbolA'] = symbolA_list
    bp_complex_i_df.loc[:,'SymbolB'] = symbolB_list
    bp_complex_i_df.loc[:,'pW'] = pW_list
    bp_complex_i_df.loc[:,'pNI'] = pNI_list
    bp_complex_i_df.loc[:,'pInt'] = pInt_list
    
    return bp_complex_i_df
    
def get_prop_edges_in_complex_identified(bp_PPI_G, Corum_DF, Complex_ID):
    '''
    Calculates proportion of all possible edges identified from BioPlex (AP-MS) 
    PPIs for a CORUM complex.
    
    This function returns the proportion of all possible PPIs identified 
    through AP-MS between the proteins in a specified CORUM complex.

    Parameters
    ----------
    DataFrame of PPIs : Pandas DataFrame
    DataFrame of CORUM complexes : Pandas DataFrame
    Corum Complex ID: int

    Returns
    -------
    Float
        The proportion of interactions between all proteins in CORUM complex 
        identified through AP-MS PPI data

    Examples
    --------
    # (1) Obtain the latest version of the 293T PPI network
    # (2) Obtain NetworkX graph representation of 293T PPI network
    # (3) Obtain CORUM complexes
    # (4) Get proportion of interactions identified for a specified CORUM complex 
    # using PPI data

    >>> bp_293t_df = getBioPlex('293T', '3.0')
    >>> bp_293t_G = bioplex2graph(bp_293t_df)
    >>> Corum_DF = getCorum()
    >>> get_prop_edges_in_complex_identified(bp_293t_G, Corum_DF, 2851)
    '''
    # store gene UNIPROT IDs that belong to this complex in a list
    genes_in_complex_i = (Corum_DF[Corum_DF.complex_id == Complex_ID].loc[:,
                                'subunits_uniprot_id'].values[0].split(';'))
    
    # get subgraph induced by the subset of nodes in this CORUM complex
    bp_complex_i_G = bp_PPI_G.subgraph(genes_in_complex_i)
        
    # create a complete graph from the nodes of complex graph 
    # (all possible interactions between proteins)
    bp_complex_i_G_complete = nx.Graph()
    bp_complex_i_G_complete.add_nodes_from(genes_in_complex_i)
    bp_complex_i_G_complete.add_edges_from(
                    itertools.combinations(genes_in_complex_i, 2))
    
    # calculate proportion of interactions between proteins in complex 
    # identified through AP-MS
    prop_edges_identified = (float(len(list(bp_complex_i_G.edges)))/
                             float(len(list(bp_complex_i_G_complete.edges))))
    
    # return proportion of edges ID'd through AP-MS, round to 3 decimal places
    return round(prop_edges_identified, 3)

def resampling_test_for_uniprot_list(bp_PPI_G, uniprot_list, 
                    num_resamples = 1000, preserve_node_degree = False):
    '''
    Run resampling test to check for enrichment of BioPlex (AP-MS) PPIs 
    for a given list of proteins (uniprot IDs).

    This function returns a p-value after running a resampling test by 
    1. taking the number of proteins in the specified list of uniprot IDs (N), 
    2. choosing N random proteins from the Graph generated by all of the 
       PPI data (G),
    3. calculating the number of edges in the Subgraph (S) induced by N random 
       proteins (with the same proportion of baits (+/- 10%) as the complex) 
       and storing this value (E_i),
    4. if preserve_node_degree option is invoked, then baits & preys in S 
       must have the same degree distribution as baits & preys in the
       network generated by given uniprot IDs, respectively,
    5. repeating steps 1-3 num_resamples times to create a null distribution, 
    6. calculating the number of edges between N proteins in the 
       complex (E), 
    7. returning a p-value by calculating the proportion of values 
       [E_1, E_2, ... , E_num_resamples] that are greater than or equal to E.

    Parameters
    ----------
    Network of PPIs : NetworkX graph
    list of uniprot IDs : list
    Number of resamplings: int
    option to preserve degree distribution in subgraphs : bool

    Returns
    -------
    Float
        A p-value from a resampling test to check for enrichment of PPIs 
        detected between proteins in list

    Examples
    --------
    # (1) Obtain the latest version of the 293T PPI network
    >>> bp_293t_df = getBioPlex('293T', '3.0')
    # (2) Obtain NetworkX graph representation of 293T PPI network
    >>> bp_293t_G = bioplex2graph(bp_293t_df)
    # (3) Obtain CORUM complexes
    >>> Corum_DF = getCorum()
    # (4) Get list of uniprots for Arp2/3 complex
    >>> UniProts_Arp_2_3 = get_UniProts_from_CORUM(Corum_DF, Complex_ID = 27)
    # (5) Calculate p-value to check for enrichment of edges in 
    #     Arp2/3 protein complex
    >>> resampling_test_for_CORUM_complex(bp_293t_G, UniProts_Arp_2_3, 1000)
    '''
    # get subgraph induced by the subset of nodes in this protein list
    bp_uniprots_i_G = bp_PPI_G.subgraph(uniprot_list)

    # number of edges detected between proteins in this CORUM 
    # complex among PPI data
    num_edges_identified_uniprot_list = float(len(list(bp_uniprots_i_G.edges)))

    # if no edges detected in this uniprot list, return an ERROR message
    if num_edges_identified_uniprot_list == 0.0:
        print('ERROR: no edges detected in PPI data for this protein list, '
              'p-value could not be computed.')

    # if at least 1 edge in complex, estimate p-value 
    # using resampling test
    else:

        # number of genes in complex (N genes)
        num_genes_in_uniprot_list = len(list(bp_uniprots_i_G.nodes))    

        # list of nodes in large network generated from PPI data
        nodes_in_overall_PPI_network = list(bp_PPI_G.nodes)
        
        # bait degrees in large PPI network
        ########################################
        # create a filter for uniprot IDs in large network that are baits
        G_baits_filter = np.array(
            [bp_PPI_G.nodes[node_i]['bait'] for node_i in bp_PPI_G.nodes])

        # get array of baits in large network
        G_baits = np.array(nodes_in_overall_PPI_network)[G_baits_filter]

        # get degree of each bait
        G_baits_degrees = bp_PPI_G.degree(G_baits)

        # prey degrees in large PPI network
        ########################################
        # create a filter for uniprot IDs in complex that are preys
        G_preys_filter = np.array(
            [not bp_PPI_G.nodes[node_i]['bait'] for node_i in bp_PPI_G.nodes])

        # get array of preys in large network
        G_preys = np.array(nodes_in_overall_PPI_network)[G_preys_filter]

        # get degree of each prey
        G_preys_degrees = bp_PPI_G.degree(G_preys)

        # bait degrees in complex
        ########################################
        # get array of uniprot IDs in complex
        bp_uniprots_i_G_nodes = np.array(bp_uniprots_i_G.nodes)

        # create a filter for uniprot IDs in complex that are baits
        bp_uniprots_i_G_baits_filter = np.array(
            [bp_uniprots_i_G.nodes[node_i]['bait'] for 
                 node_i in bp_uniprots_i_G.nodes])

        # get list of baits in complex
        bp_uniprots_i_G_baits = bp_uniprots_i_G_nodes[
            bp_uniprots_i_G_baits_filter]

        # get degree of each bait
        bp_uniprots_i_G_baits_degrees = bp_uniprots_i_G.degree(
                                        bp_uniprots_i_G_baits)

        # get degree distribution of baits
        bp_uniprots_i_G_baits_degree_distr = Counter(
            [bp_uniprots_i_G_baits_degrees[uniprot_i] for 
             uniprot_i in bp_uniprots_i_G_baits])
        
        # remove count for any "0" degree baits
        # would not contribute to edges in S
        if 0 in bp_uniprots_i_G_baits_degree_distr.keys():
            del bp_uniprots_i_G_baits_degree_distr[0]
        
        # find number of baits in complex
        
        # if preserve node degree TRUE
        # exclude baits w/ degree 0 from bait count
        # b/c random subgraphs won't have these
        if preserve_node_degree == True:
            bp_complex_i_baits_num = (np.sum(list(
                    bp_uniprots_i_G_baits_degree_distr.values())))
        
        # else count all baits in complex
        else:
            bp_complex_i_baits_num = (np.sum(bp_uniprots_i_G_baits_filter))
        
        # find proportion of baits in complex
        bp_complex_i_baits_prop = (float(bp_complex_i_baits_num) / 
                                   float(num_genes_in_uniprot_list))

        # prey degrees in complex
        ########################################
        # get array of uniprot IDs in complex
        bp_uniprots_i_G_nodes = np.array(bp_uniprots_i_G.nodes)

        # create a filter for uniprot IDs in complex that are preys
        bp_uniprots_i_G_preys_filter = np.array(
            [not bp_uniprots_i_G.nodes[node_i]['bait'] for 
             node_i in bp_uniprots_i_G.nodes])

        # get list of prey in complex
        bp_uniprots_i_G_preys = bp_uniprots_i_G_nodes[
            bp_uniprots_i_G_preys_filter]

        # get degree of each prey
        bp_uniprots_i_G_preys_degrees = bp_uniprots_i_G.degree(
            bp_uniprots_i_G_preys)

        # get degree distribution of preys
        bp_uniprots_i_G_preys_degree_distr = Counter(
            [bp_uniprots_i_G_preys_degrees[uniprot_i] for 
             uniprot_i in bp_uniprots_i_G_preys])
        
        # remove count for any "0" degree preys
        # would not contribute to edges in S
        if 0 in bp_uniprots_i_G_preys_degree_distr.keys():
            del bp_uniprots_i_G_preys_degree_distr[0]

        # list that will store number of edges detected in each subgraph
        num_edges_random_subgraphs = []

        # iterate through num_resamples random subgraphs induced by N nodes
        S_i = 0
        while S_i < num_resamples:
            
            # if degree distribution option is invoked
            if preserve_node_degree == True:

                # baits in S
                ########################################
                # preserve degrees by randomly pulling baits
                # with same degree for each uniprot in complex
                S_rando_baits = []
                for deg_i, uniprot_count_i in zip(
                    bp_uniprots_i_G_baits_degree_distr.keys(), 
                    bp_uniprots_i_G_baits_degree_distr.values()):

                    # get baits in G with same degree
                    G_baits_with_deg_i = list(
                        G_baits[np.array(
                            [G_baits_degrees[bait_i] == deg_i for 
                             bait_i in G_baits])])
                    
                    # choose N baits w/ same degree at random w/o replacement
                    N_rando_baits_from_PPI_network = random.sample(
                        G_baits_with_deg_i, uniprot_count_i)

                    # add to list of random baits
                    S_rando_baits = (S_rando_baits + 
                                     N_rando_baits_from_PPI_network)

                # preys in S
                ########################################
                # preserve degrees by randomly pulling preys
                # with same degree for each uniprot in complex
                S_rando_preys = []
                for deg_i, uniprot_count_i in zip(
                    bp_uniprots_i_G_preys_degree_distr.keys(), 
                    bp_uniprots_i_G_preys_degree_distr.values()):

                    # get preys in G with same degree
                    G_preys_with_deg_i = list(
                        G_preys[np.array([G_preys_degrees[prey_i] == deg_i for 
                                          prey_i in G_preys])])

                    # choose N preys w/ same degree at random w/o replacement
                    N_rando_preys_from_PPI_network = random.sample(
                        G_preys_with_deg_i, uniprot_count_i)

                    # add to list of random preys
                    S_rando_preys = (S_rando_preys + 
                                     N_rando_preys_from_PPI_network)

                # combine random Preys & Baits while preserving degrees
                N_rando_nodes_from_PPI_network = S_rando_baits + S_rando_preys

            # if degree distribution not invoked
            elif preserve_node_degree == False:

                # choose N genes at random without replacement
                N_rando_nodes_from_PPI_network = random.sample(
                    nodes_in_overall_PPI_network, num_genes_in_uniprot_list)

            # get subgraph induced by random subset of nodes
            bp_PPI_S = bp_PPI_G.subgraph(N_rando_nodes_from_PPI_network)

            # check to see if nodes in subgraph have the same proportion of 
            # baits as the subgraph induced by the CORUM complex
            # excluding baits that had degree 0
            bp_S_baits_num = np.sum(
                [bp_PPI_S.nodes[node_i]['bait'] for node_i in bp_PPI_S.nodes])
            bp_S_baits_prop = (float(bp_S_baits_num) / 
                               float(num_genes_in_uniprot_list))

            # proportion of baits in CORUM complex & S are the same (+/- 10%)
            if abs(bp_complex_i_baits_prop - bp_S_baits_prop) <= 0.1:

                # calculate the number of edges detected within 
                # subgraph induced by random nodes
                num_edges_S = float(len(list(bp_PPI_S.edges)))

                # store in list that contains resamplings
                num_edges_random_subgraphs.append(num_edges_S)

                # count this as a resampling
                S_i += 1

        # convert list to numpy array
        num_edges_random_subgraphs = np.array(num_edges_random_subgraphs)

        # calculate proportion of subgraphs that had more edges than edges 
        # detected in complex (p-val from resampling test)
        p_val = (float(np.sum(num_edges_random_subgraphs 
                                >= num_edges_identified_uniprot_list) + 1.0) / 
                            (float(num_resamples) + 1.0))
    return p_val

def get_interacting_chains_from_PDB(PDB_ID_structure_i, protein_structure_dir, dist_threshold):
    '''
    Retreive chain pairs that are physically close to eachother from 
    PDB structure.
    
    This function downloads the PDB structure that is specified from the input 
    PDB ID into the input directory, then computes the pairwise distances 
    between all atoms for each pair of chains in the structure. A list of 
    chain pairs that are interacting (have at least a pair of 
    atoms < dist_threshold angstroms apart) is returned.

    Parameters
    ----------
    PDB ID: str
    directory to store PDB file: str
    distance threshold: int

    Returns
    -------
    Interacting Chains
        List of chain pairs from PDB structure that have at least 
        one pair of atoms located < distance threshold apart.

    Examples
    --------
    # (1) Obtain list of interacting chains from 6YW7 structure
    >>> interacting_chains_list = get_interacting_chains_from_PDB('6YW7', '.', 6)
    Downloading PDB structure '6YW7'...
    >>> interacting_chains_list
    [['A', 'D'], ['A', 'E'], ['A', 'B'], ['D', 'F'], ['B', 'F'], ['B', 'G'], ['F', 'G'], ['F', 'C'], ['G', 'C']]
    '''
    # download structure from PDB
    pdbl = PDBList()
    PBD_file_path = pdbl.retrieve_pdb_file(PDB_ID_structure_i, 
                                           pdir=protein_structure_dir, 
                                           file_format='pdb', 
                                           overwrite=True)

    # create a structure object
    parser = PDBParser()
    structure = parser.get_structure(PDB_ID_structure_i, PBD_file_path)

    model = structure[0]
    chain_IDs = [chain.get_id() for chain in model] # get a list of all chains

    # we want to test every pair of chains to see if they have any atoms 
    # that are < 6 angstroms in distance
    possible_chain_pairs = list(itertools.combinations(chain_IDs, 2))

    chain_pairs_direct_interaction = []
    # iterate through all chain pairs and check to see if any atoms are close
    for chain_i_id, chain_j_id in possible_chain_pairs:

        # get chain objects from models
        chain_i = model[chain_i_id]
        chain_j = model[chain_j_id]

        # get all atoms from each chain, 'A' stands for ATOM
        atom_list_i = Selection.unfold_entities(chain_i, "A")
        atom_list_j = Selection.unfold_entities(chain_j, "A")

        # get the coordinates for the atom in each chain
        atom_coords_i = ([atom_list_i[k].coord 
                          for k in range(0,len(atom_list_i))])
        atom_coords_i = np.vstack(atom_coords_i)

        atom_coords_j = ([atom_list_j[k].coord 
                          for k in range(0,len(atom_list_j))])
        atom_coords_j = np.vstack(atom_coords_j)

        # compute pairwise distances betweeen all atoms from different chains
        dists = cdist(atom_coords_i, atom_coords_j)

        # if a pair of atoms < dist_threshold angstroms apart, store as interacting chains
        if np.sum(dists < dist_threshold) >= 1:
            chain_pairs_direct_interaction.append([chain_i_id, chain_j_id])
        
    return chain_pairs_direct_interaction

def make_request(url, mode, pdb_id):
    '''
    Make requests to PDBe API.
    
    This function can make GET and POST requests to the PDBe API.
    
    Parameters
    ----------
    url: str
    mode: str
    pdb_id: str
    
    Returns
    -------
        JSON or None
    '''
    if mode == "get":
        response = requests.get(url=url+pdb_id)
    elif mode == "post":
        response = requests.post(url, data=pdb_id)

    if response.status_code == 200:
        return response.json()
    else:
        print("[No data retrieved - %s] %s" 
              % (response.status_code, response.text))
    
    return None

def get_mappings_data(pdb_id):
    '''
    Get mappings data for PDB ID.
    
    This function will retreive the mappings data from
    the PDBe API using the make_request() function.
    
    Parameters
    ----------
    pdb_id: str
    
    Returns
    -------
        JSON of mappings or None
    '''
    # specify URL
    base_url = "https://www.ebi.ac.uk/pdbe/"
    api_base = base_url + "api/"
    uniprot_mapping_url = api_base + 'mappings/uniprot/'
    
    # Check if the provided PDB id is valid
    # There is no point in making an API call
    # with bad PDB ids
    if not re.match("[0-9][A-Za-z][A-Za-z0-9]{2}", pdb_id):
        print("Invalid PDB id")
        return None
    
    # GET the mappings data
    mappings_data = make_request(uniprot_mapping_url, "get", pdb_id)
    
    # Check if there is data
    if not mappings_data:
        print("No data found")
        return None
    
    return mappings_data

def list_uniprot_pdb_mappings(pdb_id):
    '''
    Get PDB chain to UniProt mappings.
    
    This function retrieves PDB > UniProt mappings using the 
    get_mappings_data() function, the parses the resulting 
    JSON to construct a dictionary where each key is a chain
    from the PDB structure, and the corresponding value for
    each is a list of UniProt IDs that map to the chain from 
    the SIFTS project.
    
    Parameters
    ----------
    pdb_id: str
    
    Returns
    -------
    Chain to UniProt Map
        Dictionary of PDB ID chain to UniProt ID mappings
        
    Examples
    --------
    # (1) Obtain a mapping of PDB ID 6YW7 chains to UniProt IDs
    >>> chain_to_UniProt_mapping_dict = list_uniprot_pdb_mappings('6YW7')
    >>> len(chain_to_UniProt_mapping_dict)
    7
    >>> sorted(chain_to_UniProt_mapping_dict.keys())
    ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    >>> sorted(chain_to_UniProt_mapping_dict.values())
    [['O15144'], ['O15145'], ['O15511'], ['P59998'], ['P61158'], ['P61160'], ['Q92747']]
    '''
    # convert to PDB id to lower case
    pdb_id = pdb_id.lower()
    
    # Getting the mappings data
    mappings_data = get_mappings_data(pdb_id)
    
    # If there is no data, return None
    if not mappings_data:
        return None
    
    # dictionary that stores UniProt > chain id mappings
    uniprot_chain_mapping_dict = {}
    
    uniprot = mappings_data[pdb_id]["UniProt"]
    for uniprot_id in uniprot.keys():
        mappings = uniprot[uniprot_id]["mappings"]
        
        # store the chain ids that correspond to this UniProt ID
        uniprot_chain_mapping_dict[uniprot_id] = []
        
        for mapping in mappings:
            entity_id = mapping["entity_id"]
            
            chain_id = mapping["chain_id"]
            uniprot_chain_mapping_dict[uniprot_id].append(chain_id)
            
            pdb_start = mapping["start"]["residue_number"]
            pdb_end = mapping["end"]["residue_number"]
            uniprot_start = mapping["unp_start"]
            uniprot_end = mapping["unp_end"]
        
    # "flip" the uniprot > pdb chain mapping
    # get all unique chain IDs
    chain_IDs = (list(set([item for sublist in 
                           uniprot_chain_mapping_dict.values() for item 
                           in sublist])))

    chain_uniprot_mapping_dict = {}
    # iterate through every chain ID
    for chain_i in chain_IDs:

        # iterate through every uniprot ID and check to 
        # see if chain ID is mapped
        chain_uniprot_mapping_dict[chain_i] = []
        for uniprot_i in uniprot_chain_mapping_dict.keys():
            if chain_i in uniprot_chain_mapping_dict[uniprot_i]:
                chain_uniprot_mapping_dict[chain_i].append(uniprot_i)
                
    return chain_uniprot_mapping_dict

def PDB_chains_to_uniprot(interacting_chains_list, 
                          chain_to_UniProt_mapping_dict):
    '''
    Get interacting chains from PDB structure mapped to UniProt IDs.
    
    This function takes the list of interacting chains from function
    get_interacting_chains_from_PDB() and the chain to UniProt mappings
    from function list_uniprot_pdb_mappings() and returns a list of 
    interacting chains using UniProt IDs.
    
    Parameters
    ----------
    Interacting Chains: list
    Chain to UniProt Map: dict
    
    Returns
    -------
    Interacting Chains
        List of interacting chains using UniProt IDs
    
    Examples
    --------
    # (1) Obtain list of interacting chains from 6YW7 structure
    # (2) Obtain a mapping of PDB ID 6YW7 chains to UniProt IDs
    # (3) Obtain list of interacting chains from 6YW7 
    #     structure using UniProt IDs 
    >>> interacting_chains_list = get_interacting_chains_from_PDB('6YW7', '.', 6)
    Downloading PDB structure '6YW7'...
    >>> chain_to_UniProt_mapping_dict = list_uniprot_pdb_mappings('6YW7') 
    >>> interacting_UniProt_IDs = PDB_chains_to_uniprot(interacting_chains_list, chain_to_UniProt_mapping_dict)
    >>> sorted(interacting_UniProt_IDs)
    [['O15144', 'P59998'], ['O15511', 'Q92747'], ['P59998', 'O15511'], ['P59998', 'Q92747'], ['P61158', 'O15144'], ['P61158', 'O15145'], ['P61158', 'P61160'], ['P61160', 'O15511'], ['P61160', 'P59998']]
    '''
    interacting_UniProt_IDs = []
    for interacting_chain_pair_i in interacting_chains_list:

        # get UniProt IDs that map to each chain ID
        chain_i_UniProts = (
            chain_to_UniProt_mapping_dict[interacting_chain_pair_i[0]])
        chain_j_UniProts = (
            chain_to_UniProt_mapping_dict[interacting_chain_pair_i[1]])

        # store every pair of UniProt IDs that correspond 
        # to the interacting chains
        for chain_i_UniProt_ID in chain_i_UniProts:
            for chain_j_UniProt_ID in chain_j_UniProts:
                interacting_UniProt_IDs.append(
                    [chain_i_UniProt_ID,chain_j_UniProt_ID])
                
    return interacting_UniProt_IDs

def PDB_to_interacting_chains_uniprot_maps(PDB_ID, 
                                           protein_structure_dir, 
                                           interact_dist_threshold):
    '''
    Get interacting chains from PDB structure mapped to UniProt IDs and 
    PDB chain to UniProt mappings.
    
    This is a wrapper function for functions 
    (1) get_interacting_chains_from_PDB(),
    (2) list_uniprot_pdb_mappings(), and
    (3) PDB_chains_to_uniprot() 
    to get a list of interacting chains from PDB structure using UniProt labels 
    and the chain-to-UniProt mapping for this PDB structure.
    
    Parameters
    ----------
    PDB ID: str
    directory to store PDB file: str
    distance threshold: int
    
    Returns
    -------
    Chain to UniProt Map
        Dictionary of PDB ID chain to UniProt ID mappings
    Interacting Chains
        List of interacting chains using UniProt IDs
    
    Examples
    --------
    >>> (chain_to_UniProt_mapping_dict, interacting_UniProt_IDs) = PDB_to_interacting_chains_uniprot_maps('6NMI', '.', 6)
    Downloading PDB structure '6NMI'...
    '''
    # get list of chain pairs that interact in PDB structure 
    interacting_chains_list = get_interacting_chains_from_PDB(PDB_ID, 
                                                    protein_structure_dir, 
                                                    interact_dist_threshold)
    
    # get chain > UniProt ID mappings for this PDB structure
    chain_to_UniProt_mapping_dict = list_uniprot_pdb_mappings(PDB_ID)
    
    # get list of chains pairs that interact in PDB structure using UniProts
    interacting_UniProt_IDs = PDB_chains_to_uniprot(interacting_chains_list, 
                                                chain_to_UniProt_mapping_dict)

    return [chain_to_UniProt_mapping_dict, interacting_UniProt_IDs]