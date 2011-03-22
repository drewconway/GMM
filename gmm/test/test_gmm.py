#!/usr/bin/env python
# encoding: utf-8
"""
test_gmm.py

Purpose:  Tests for base gmm class

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2011-03-22

"""

import unittest
import networkx as nx
import gmm

class test_gmm(unittest.TestCase):
    """Tests for base gmm class"""
            
    # GMM object place holder
    full_model=None     
    
    # Edges of a five-cycle graph
    cycle_edges=[(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]
    
    # Base graph
    five_cycle=nx.Graph(data=cycle_edges)               
    
    # Test graph
    test_triangle=nx.Graph(data=[(0,1),(1,2),(0,2)])
    
    def setUp(self):
        """Tests that all models from setup pass inspection."""
        
        # Node ceiling termination rule for 100 nodes
        def node_ceiling(G):          
            if G.number_of_nodes()>=100:
                return True 
            else:
                return False

        # Simple random growth rule: connects random node from base to random node from new
        def rand_add(base, new):
            from numpy.random import randint
            new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
            new_base=nx.compose(base,new)
            base_connector=randint(base.number_of_nodes())
            new_connector=randint(min(new.nodes()),max(new.nodes()))
            new_base.add_edge(base_connector,new_connector)
            return new_base
        
        # This models should pass without exception
        self.full_model=gmm.gmm(self.five_cycle,T=node_ceiling,R=rand_add)
        
    def testBase(self):
        """Test that gmm objects can have base graphs altered"""
        # First test that set and get functions work
        self.assertEquals(self.full_model.get_base().edges(),self.cycle_edges)
        new_edges=[(0,1),(0,2),(1,2)]
        new_edges.sort()
        self.full_model.set_base(nx.Graph(data=new_edges))
        self.assertEquals(self.full_model.get_base().edges(),new_edges)
        # Test that reversion works
        self.full_model.revert_base()
        self.assertEquals(self.full_model.get_base().edges(),self.cycle_edges)
    
    def testTermination(self):
        """Test that termination rule is set correctly, and works"""
        self.full_model.set_base(nx.complete_graph(100))
        self.assertTrue(self.full_model.apply_termination())
        self.full_model.revert_base()
        self.assertFalse(self.full_model.apply_termination())
    
    def testRule(self):
        """Test that growth rule is set correctly, and works"""
        # First test that it worked as expected
        self.assertEquals(self.full_model.apply_rule(self.test_triangle).number_of_nodes(),8)
        # Test that it is setting base correctly when set_result is True
        self.full_model.apply_rule(self.test_triangle,set_result=True)
        new_nodes=self.full_model.get_base().nodes()
        new_nodes.sort()
        self.assertEquals(new_nodes,range(8))
        self.full_model.revert_base()
        # Test that coercion is working when new structure type does not match base
        directed_triangle=self.test_triangle.to_directed()
        self.full_model.apply_rule(directed_triangle,set_result=True)
        self.assertFalse(self.full_model.get_base().is_directed())
    
if __name__ == '__main__':
    unittest.main()