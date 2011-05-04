#!/usr/bin/env python
# encoding: utf-8
"""
erdos_renyi.py

Purpose:  Recover the Erdos-Renyi binomial random graph model using GMM
          see, Erdős and A. Rényi, On Random Graphs, Publ. Math. 6, 290 (1959).

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2011-03-23

Copyright (c) 2011, under the Simplified BSD License.  
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import sys
import os
import networkx as nx
import gmm
from numpy import arange, random
from scipy import stats
import matplotlib.pylab as plt

def binomial_growth(base, new):
    """
    For each node in new_nodes, add edge to nodes in base_nodes with probability p. 
    """
    p=0.5
    
    # To keep new nodes from over-writing current ones rename the new nodes starting
    # from the last node in base
    new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
    new_nodes=new.nodes()
    base_nodes=base.nodes()
    new_base=nx.compose(base,new)

    # Add new structure to base graph given a binomial probability of 
    # forming a tie between a  node from the new structure and 
    # each node in the base structure
    for n in new_nodes:
        edge_test=zip(random.uniform(size=len(base_nodes)), base_nodes)
        for d,m in edge_test:
            if (d <= p):
                new_base.add_edge(m, n)
    return new_base
    
def binomial_simulation(graph_set, seed=None, verbose=False):
    """
    Given a set of binomial random graphs, simulate using above growth rule
    """
    
    simulated_graphs=list()
    
    # Iterate over graph set to produce several simulations
    for i in xrange(0,len(graph_set)-1):
        # Required to change the growth rule dynamically, 
        # so we define it inside the simulation's for-loop
        termination_size=graph_set[i+1].number_of_nodes()
        
        def node_ceiling(G):
            if G.number_of_nodes()>=termination_size:
                return False
            else:
                return True
        
        # Setup GMM 
        erdos_renyi_model=gmm.gmm(graph_set[i])
        erdos_renyi_model.set_termination(node_ceiling)
        erdos_renyi_model.set_rule(binomial_growth)
        sim_name="GMM SIMULATION-"+erdos_renyi_model.get_base().name
        gmm.algorithms.simulate(erdos_renyi_model,tau=3, poisson=True, seed=seed, new_name=sim_name)
        
        # Run simulation
        simualted_erdos_renyi=erdos_renyi_model.get_base()
        
        # Print the results to stdout
        if verbose:
            print(nx.info(graph_set[i+1]))
            print(nx.info(simualted_erdos_renyi))
            print("")
        simulated_graphs.append(simualted_erdos_renyi)
    
    # Return a list of simulated graphs
    return(simulated_graphs)

    
def main():
    
    """
    1. Generate a set of Erdos-Renyi graphs to experiment on
    """
    
    # Set seed for random graph
    random_seed=851982
    
    # Current working directory
    cwd=os.getcwd()
    
    # Produce a set of Erdos-Renyi binomial random graphs of sizes ranging from 25 to 100 nodes.
    # For the purposes of this exerpiment the value p=0.5 is constant.  As such, it is hard-coded 
    # in simulation process for both the base graphs and the GMM simulations.
    graph_set=map(lambda i: nx.erdos_renyi_graph(i, p=0.5, seed=random_seed), arange(25,125,25))
    
    # Output the base graphs as GraphML
    for g in graph_set:
        nx.write_graphml(g, cwd+"/erdos_renyi/data/"+g.name+"_"+str(g.number_of_nodes())+".graphml") 
    
    print("Base Erdos-Renyi graphs generated")
    
    """
    2. Give set of graphs, simulate graph_set[i] to graph_set[i+1] with GMM using random graph rule
    and node ceiling termination given graph_set[i+1].number_of_nodes()
    """
    gmm_simulations=binomial_simulation(graph_set, seed=random_seed, verbose=True)
    
    # Output the simulations as GraphML 
    for g in gmm_simulations:
        nx.write_graphml(g, cwd+"/erdos_renyi/data/"+g.name+"_"+str(g.number_of_nodes())+".graphml")

if __name__ == '__main__':
	main()

