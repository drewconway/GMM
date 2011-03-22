#!/usr/bin/env python
# encoding: utf-8
"""
zachary_regen.py

Purpose:  Attempt to regenerate very Zachary's Karate Club data (very small). Will begin with full data 
            (34 nodes), do some random deletion, and then attempt to regenerate with simple growth rule 
            and node uses node ceiling as termiantion rule.
            
            Zachary network data available from Mark Newman at:
            http://networkdata.ics.uci.edu/data.php?id=105
            
            Original citiation:
            W. W. Zachary, "An information flow model for conflict and fission in small groups", 
            Journal of Anthropological Research 33, 452-473 (1977).

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2010-11-04

Copyright (c) 2010, under the Simplified BSD License.  
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import sys
import os
import matplotlib.pylab as plt
import networkx as nx
import gmm
from numpy.random import random_integers, uniform

# Randomly delete some number of nodes from a given graph
def rand_delete(G, num_nodes):
    G=nx.convert_node_labels_to_integers(G,first_label=0)
    nodes_to_delete=list(random_integers(low=0,high=len(G.nodes()),size=num_nodes))
    G.remove_nodes_from(nodes_to_delete)
    isos=nx.isolates(G)
    G.remove_nodes_from(isos)
    return(G)

# Termnation rule: Attempt to match orignal 34 nodes
def node_ceiling_34(G):
    if G.number_of_nodes()>=34:
        return False
    else:
        return True
        
# Termnation rule: grow original network to have at least 100 members
def node_ceiling_100(G):
    if G.number_of_nodes()>=100:
        return False
    else:
        return True        
        
# Growth rule for karate club
def karate_rule(base, new):
    """
    The original karate club has two densly connected clusters, with a few critical bridges. This rule attempts 
    to simulate a growth function from this by first bringing together all components, then with a simple
    preferential attachment model to those actors with relatively high Eigenvector centrality.  The original 
    Zachary data forms two tightly clustered communities with a select few bridges. This growth rule is meant to 
    re-constrcut this in reverse, by first forming the bridges, then the polar hubs.
    """
    # If the base graph has multiple components, first attempt to unify them into a single component
    if nx.components.number_connected_components(base)>1:
        comps=nx.connected_components(base)
        # Select two random components and form bridge with new structure
        rand_comps=random_integers(low=0,high=len(comps)-1,size=2)
        # Select random node from each component
        rand0=comps[rand_comps[0]][random_integers(low=0,high=len(comps[rand_comps[0]])-1)]
        rand1=comps[rand_comps[1]][random_integers(low=0,high=len(comps[rand_comps[1]])-1)]
        while rand0==rand1:
            rand1=comps[rand_comps[1]][random_integers(low=0,high=len(comps[rand_comps[1]])-1)]
        outer_bound=[rand0,rand1]
        outer_bound.extend(range(base.number_of_nodes(),base.number_of_nodes()+((new.number_of_nodes())-1)))
        mapping=dict(zip(new.nodes(),outer_bound))
        new=nx.relabel_nodes(new,mapping)
    else:
        # Use Eigenvector centrality as pref attachment attribute
        cent=nx.eigenvector_centrality_numpy(base).items()
        # Normalize values to sum to 1.0
        norm_const=sum([(b) for (a,b) in cent])
        pref_prob=[(b/norm_const) for (a,b) in cent]
        # Step through probability mass to find a node to attach to. Same method used in 
        # gmm.algorithms.draw_structure to select a probability weighted motif from the set.
        draw=uniform()
        node_index=0
        mass_sum=pref_prob[node_index]
        while draw>mass_sum:
            node_index+=1
            mass_sum+=pref_prob[node_index]
        rand0=cent[node_index][0] # Return the appropriate node ID
        outer_bound=[rand0]
        outer_bound.extend(range(base.number_of_nodes(),base.number_of_nodes()+((new.number_of_nodes())-1)))
        mapping=dict(zip(new.nodes(),outer_bound))
        new=nx.relabel_nodes(new,mapping)
    return nx.compose(base,new)
    


def main():
    # Load Zachary data, randomly delete nodes, and report
    zachary=nx.Graph(nx.read_pajek("karate.net")) # Do not want graph in default MultiGraph format
    zachary.name="Original Zachary Data"
    print(nx.info(zachary))
    zachary_subset=rand_delete(zachary, 15) # Remove half of the structure
    zachary_subset.name="Randomly Deleted Zachary Data"
    print(nx.info(zachary_subset))
    
    # Create model, and simulate
    zachary_model=gmm.gmm(zachary_subset,R=karate_rule,T=node_ceiling_34)
    gmm.algorithms.simulate(zachary_model,4,poisson=False,new_name="Simulation from sample")  # Use tau=4 because data is so small (it's fun!)
    
    # Report and visualize
    print(nx.info(zachary_model.get_base()))
    fig=plt.figure(figsize=(30,10))
    fig.add_subplot(131)
    nx.draw_spring(zachary,with_labels=False,node_size=45,iterations=5000)
    plt.text(0.01,-0.1,"Original Karate Club",color="darkblue",size=20)
    fig.add_subplot(132)
    nx.draw_spring(zachary_subset,with_labels=False,node_size=45,iterations=5000)
    plt.text(0.01,-0.1,"Random sample of Karate Club",color="darkblue",size=20)
    fig.add_subplot(133)
    nx.draw_spring(zachary_model.get_base(),with_labels=False,node_size=45,iterations=5000)
    plt.text(0.01,-0.1,"Simulation from random sample",color="darkblue",size=20)
    plt.savefig("zachary_simulation.png")


if __name__ == '__main__':
	main()

