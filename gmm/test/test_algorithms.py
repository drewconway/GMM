#!/usr/bin/env python
# encoding: utf-8
"""
test_algorithms.py

Purpose:  Test for all algorithms used to simulate network structure from
          graph motifs.

Author:   Drew Conway
Email:    drew.conway@nyu.edu

"""

import unittest
import copy
import networkx as nx
import gmm

class test_algorithms(unittest.TestCase):
    
    # GMM object place holder
    base_model=None 
    base_directed=None    

    # Base graph
    five_cycle=nx.cycle_graph(5)

    # Triangle for motif checks
    test_triangle=nx.Graph(data=[(0,1),(1,2),(0,2)])
    directed_triangle=copy.deepcopy(test_triangle).to_directed()

    # Test tau value
    test_tau=3

    def setUp(self):
        # Most basic GMM; no growth or termination rules.
        self.base_model=gmm.gmm(self.five_cycle)
        self.base_directed=gmm.gmm(copy.deepcopy(self.five_cycle).to_directed())
        
    def test_is_gmm(self):
        """Tests that functions return appropriate logical value"""
        self.assertTrue(self.base_model.am_gmm)
    
    def test_motif_counts(self):
        """Tests that both undirected and directed motif counts are returned correctly"""
        base_counts=gmm.algorithms.motif_counts(self.base_model,self.test_tau)
        self.assertEquals(sum([(c) for (a,b,c) in base_counts]),20)
        directed_counts=gmm.algorithms.motif_counts(self.base_directed,self.test_tau)
        self.assertEquals(sum([(c) for (a,b,c) in directed_counts]),10)
    
    def test_get_motifs(self):
        """Test that the appropriate graph motifs are returned given tau"""
        base_motifs=gmm.algorithms.get_motifs(self.test_tau,False)
        directed_motifs=gmm.algorithms.get_motifs(self.test_tau,True)
        self.assertTrue(len(base_motifs)==3)
        self.assertTrue(len(directed_motifs)==8)
        
    def test_draw_structure(self):
        """Checks that the random selection of structure is implemented properly"""
        test_motifs=[(0,False,0.0),(1,True,1.0),(2,False,0.0)]
        non_random_draw=gmm.algorithms.draw_structure(test_motifs)
        self.assertTrue(non_random_draw)
        
    def test_poisson(self):
        """Tests that the Poisson probability mass is returned correctly for
        some set of counts
        """
        base_mass=gmm.algorithms.poisson_mass(gmm.algorithms.motif_counts(self.base_model,self.test_tau))
        base_positive=map(lambda x: x>0, [(c) for (a,b,c) in base_mass])
        self.assertEquals(base_positive.count(False),0)
        directed_mass=gmm.algorithms.poisson_mass(gmm.algorithms.motif_counts(self.base_directed,self.test_tau))
        directed_positive=map(lambda x: x>0, [(c) for (a,b,c) in directed_mass])
        self.assertEquals(directed_positive.count(False),0)
    
    def test_all_graphs(self):
        """Tests that function returns correct graphs"""
        all_base=gmm.algorithms.all_graphs(self.test_tau,False)
        all_directed=gmm.algorithms.all_graphs(self.test_tau,True)
        self.assertTrue(len(all_base)!=len(all_directed))
        self.assertTrue(len(all_base)==2)       # For tau=3 and undirected graphs there should be 2 single component graphs
        self.assertTrue(len(all_directed)==6)   # For tau=3 and directed graphs there should be 6 single component graphs

    
if __name__ == '__main__':
	unittest.main()