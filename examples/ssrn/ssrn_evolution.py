#!/usr/bin/env python
# encoding: utf-8
"""
ssrn_evolution.py

Purpose:  Create and test a GMM model for the evolution of the SSRN Conflict Studies
            eJournal from 2008 to 2009.

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2010-12-27

Copyright (c) 2010, under the Simplified BSD License.  
 For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import sys
import os
import matplotlib.pylab as plt

# Will need to edited to the path where you have saved the gmm folder
# Will be fixed once installation scripts are completed.
sys.path.append("/Users/agconway/Documents/GMM/")
from gmm import *


# A model of SSRN co-authorship evolution
def ssrn_evo(base, new):
    from numpy.random import randint,uniform
    new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
    # Growth rule: use new structure to fuse together components
    def comp_fuse(base, new):
        new_base=nx.compose(base,new)
        base_comps=nx.weakly_connected_component_subgraphs(base)
        # Select nodes to connect (fusion)
        rand_comps=randint(0,nx.number_weakly_connected_components(base),2)
        rand_node1=base_comps[0].nodes()[randint(0,base_comps[0].number_of_nodes())]
        rand_node2=base_comps[1].nodes()[randint(0,base_comps[1].number_of_nodes())]
        rand_node3=new.nodes()[randint(0,new.number_of_nodes())]
        rand_node4=new.nodes()[randint(0,new.number_of_nodes())]
        # Connect nodes
        new_base.add_edge(rand_node1,rand_node3)
        new_base.add_edge(rand_node2,rand_node4)
        return new_base
    # Stochastically fuse components (14% prob); otherwise,
    # simply add as new network component
    if 0.14 > uniform():
        new_base=comp_fuse(base,new)
        # Sanity check to maintain necessary bipartite structure
        while(nx.is_bipartite(new_base.to_undirected()) is False ):
            new_base=comp_fuse(base,new)
        return(new_base)
    else:
        return nx.compose(base,new)
    

# Termnation rule: Vertex ceiling for 2009 network of 1,446
def node_ceiling(G):
    if G.number_of_nodes()>=1446:
        return False
    else:
        return True

def main():
    # Import the network data
    ssrn_2008=nx.DiGraph(nx.read_graphml("data/ssrn_2008.graphml"))
    ssrn_2009=nx.DiGraph(nx.read_graphml("data/ssrn_2009.graphml"))
    
    # Create gmm model
    ssrn_model=gmm.gmm(ssrn_2008)
    ssrn_model.set_termination(node_ceiling)
    ssrn_model.set_rule(ssrn_evo)
    
    # Run model
    gmm.algorithms.simulate(ssrn_model,5,poisson=False)
    
    # View results
    sim_2009=ssrn_model.get_base()
    print "Estimated 2009"
    print(nx.info(sim_2009))
    print("Bipartite: "+str(nx.is_bipartite(sim_2009.to_undirected())))
    print ""
    print "Actual 2009"
    print(nx.info(ssrn_2009))
    
    # Output results
    nx.write_graphml(sim_2009, "data/sim_2009.graphml")


if __name__ == '__main__':
    main()

