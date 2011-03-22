#!/usr/bin/env python
# encoding: utf-8
"""
basic_model.py


Purpose:  Present a basic gmm model with results. This is completley niave and meant only
            as an exercise in creating a first gmm and exploring the results.  First we
            run a basic model using an undirected graph, then repeat with different 
            settings and a directed base graph.

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2010-11-01

Copyright (c) 2010, under the Simplified BSD License.  
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import sys
import os
import matplotlib.pylab as plt
import networkx as nx
import gmm

# Growth rule: randomly add new structure
def rand_add(base, new):
    from numpy.random import randint
    new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
    new_base=nx.compose(base,new)
    base_connector=randint(base.number_of_nodes())
    new_connector=randint(min(new.nodes()),max(new.nodes()))
    new_base.add_edge(base_connector,new_connector)
    return new_base
    
# Termnation rule: 250 vertex ceiling
def node_ceiling(G):
    if G.number_of_nodes()>=250:
        return False
    else:
        return True

def main():
    
    ### Undirected graph ###
    
    # Initialize model using the Petersen graph
    model=gmm.gmm(nx.petersen_graph())
    old_graph=model.get_base()
    model.set_termination(node_ceiling)
    model.set_rule(rand_add)
    
    # Run simualation with tau=4 and Poisson density for motifs
    gmm.algorithms.simulate(model,4)   

    # View results
    new_graph=model.get_base()
    print(nx.info(new_graph))
    
    # Draw graphs
    old_pos=nx.spring_layout(old_graph)
    new_pos=nx.spring_layout(new_graph,iterations=2000)
    fig1=plt.figure(figsize=(15,7))
    fig1.add_subplot(121)
    #fig1.text(0.1,0.9,"Base Graph")
    nx.draw(old_graph,pos=old_pos,node_size=25,with_labels=False)
    fig1.add_subplot(122)
    #fig1.text(0.1,0.45,"Simulation Results")
    nx.draw(new_graph,pos=new_pos,node_size=20,with_labels=False)
    fig1.savefig("undirected_model.png")
    
    ### Directed graph ###
    
    # Initialize model using random directed Barabasi-Albert model
    directed_base=nx.barabasi_albert_graph(25,2).to_directed()
    directed_model=gmm.gmm(directed_base)
    directed_model.set_termination(node_ceiling)
    directed_model.set_rule(rand_add)
    
    # Run simualation with tau=4 and Poisson density for motifs
    gmm.algorithms.simulate(directed_model,4)
    
    # View results
    new_directed=directed_model.get_base()
    print(nx.info(new_directed))
    
    # Draw directed graphs
    old_dir_pos=new_pos=nx.spring_layout(directed_base)
    new_dir_pos=new_pos=nx.spring_layout(new_directed,iterations=2000)
    fig2=plt.figure(figsize=(7,10))
    fig2.add_subplot(211)
    fig2.text(0.1,0.9,"Base Directed Graph")
    nx.draw(directed_base,pos=old_dir_pos,node_size=25,with_labels=False)
    fig2.add_subplot(212)
    fig2.text(0.1,0.45, "Simualtion Results")
    nx.draw(new_directed,pos=new_dir_pos,node_size=20,with_labels=False)
    fig2.savefig("directed_model.png")
    
    # Export files
    nx.write_graphml(model.get_base(), "base_model.graphml")
    nx.write_graphml(directed_model.get_base(), "directed_model.graphml")
    nx.write_graphml(nx.petersen_graph(), "petersen_graph.graphml")
    
    

if __name__ == '__main__':
	main()

