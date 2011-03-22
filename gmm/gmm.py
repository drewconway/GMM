#!/usr/bin/env python
# encoding: utf-8
"""
gmm.py

Purpose:  Classes used to generate networks from the graph motif model proposed
          in paper entitled, "Graph Motif Models." This method is an alternate
          process for generating network structure from observed structure.

Author:   Drew Conway
Email:    drew.conway@nyu.edu

"""

__author__="Drew Conway (drew.conway@nyu.edu)"
__docformat__ = "restructuredtext en"

import copy
import networkx as nx

class gmm(object):
    """
    This is the base class used to simulate networks with graph motif models (GMM). The 
    class inherits graph structure from NetworkX (http://networkx.lanl.gov/), a Python 
    library for the manipulation of graph data. This class contains only base data
    necessary calls to to manipulate gmm.
    
    Algorithms to simulate network structre from this object are contained in 
    gmm_simulate.py.
    
    Parameters
    ----------
    G : base graph, NetworkX Graph or DiGraph object, required at initialization
        The base graph must be a NetworkX Graph or DiGraph object with more than 
        a singleedge. It is used as the initial structure for the GMM, and is the 
        only parameter required to initialze a gmm object.
        
        **NOTE**: Nodes in base graph are not required to be integers, but if the 
        are not they are converted to integers, and the original labels are stored 
        as a dictionary by nx.convert_nodel_labels_to_integers
        
    T : model termination rule, function, optional at initialization
        The rule by which the model will terminate. Must be a function that can 
        operate on a NetworkX Graph or DiGraph object, and takes a single graph 
        object as its only parameter.
        
    R : model growth rule, function, optional at initialization
        The rule by which new structure is added to the base graph. Must be a 
        function that can operate on NetworkX Graph or DiGraph objects.  Function 
        must take exactly two arguments: 1) gmm base structure; 2) new structure 
        to be added, as NetworkX Graph or DiGraph object.  If second argument does 
        not match first it will be coerced to match.
            
    Notes
    ------
    Any deviation in function design for T or R will cause modeling errors. Basic 
    checks are performed by the gmm class at initialization, but design bugs in 
    either F or T may go undetected despite successful gmm initialization.
        
    Examples
    ----------
    # Create most basic GMM object with five node cycle graph as base.
    
    >>> import gmm
    >>> G=nx.cycle_graph(5)
    >>> model=gmm(G)
    
    # Using five node cycle graph, create gmm object with node ceiling
    # termination rule of 100 nodes
    
    >>> def degree_ceiling(G):
       ...:     if G.number_of_nodes()>=100:
       ...:         return True
       ...:     else:
       ...:         return False
    >>> model=gmm.gmm(G,degree_ceiling)
    
    # Next, add a simple random growth rule
    
    >>> def rand_add(base, new):
       ...:     from numpy.random import randint
       ...:     new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
       ...:     new_base=nx.compose(base,new)
       ...:     base_connector=randint(base.number_of_nodes())
       ...:     new_connector=randint(min(new.nodes()),max(new.nodes()))
       ...:     new_base.add_edge(base_connector,new_connector)
       ...:     return new_base
    >>> model.set_rule(rand_add)
    
    # Run simualation with tau=4 and Poisson density for motifs
    
    >>> gmm.algorithms.simulate(model,4)
    
    # View results
    
    >>> new_graph=model.get_base()
    >>> print(nx.info(new_graph))
    
    """
    ### Initialize gmm object
    def __init__(self, G, T=None,R=None):
        # Degenerate graph for testing gmm growth rule   
        self.test_graph=nx.Graph(data=[(0,1),(1,2)])           # Dyad
        self.am_gmm=True
        # Initialize GMM with base structure
        if type(G)==type(nx.Graph()) or type(G)==type(nx.DiGraph()):
            if(G.number_of_edges()>1):
                G=nx.convert_node_labels_to_integers(G,discard_old_labels=False)
                self.base=G
                self.original=copy.deepcopy(G)  # copy of graph to remain unaltered by simulations
            else:
                raise ValueError("Base graph must have at least two edges")
        else:
            raise TypeError("Base graph to gmm must be a NetworkX Graph or DiGraph object.")
        # Store termination rule if passed by user, test that it is compatible with base graph
        if T is not None:
            try:
                # Test if T is a NetworkX compatible function 
                result=T(self.base)
                if(result is True or result is False):
                    self.termination=T
                else:
                    raise TypeError("Termination rule must return boolean")
            except TypeError:
                print("T not compatible with base graph, rule set to None")
                self.termination=None
        else:
            self.termination=None
        # Store optimization rule if passed by user, test that it is compatible with base and takes 
        # NetworkX graph as second argument
        if R is not None:
            try:
                R(self.base,self.test_graph)
                self.rule=R
            except TypeError:
                print("R must be a function compatible with NetworkX graph objects, growth rule set to None.")
                self.rule=None
        else:
            self.rule=None
            
    ### gmm support functions    
    def get_base(self, original=False):
        """Returns base graph of gmm
        
        If original is True, return copy of base graph passed at initialization.  If original is False, 
        return current base graph.
        """
        if original is True:
            return self.original
        else:
            return self.base
        
    def set_base(self, G):
        """Set new base graph for gmm, but does not alter original copy"""
        if type(G)==type(nx.Graph()) or type(G)==type(nx.DiGraph()):
            if(G.number_of_edges()>1):
                G=nx.convert_node_labels_to_integers(G,discard_old_labels=False)
                self.base=G
            else:
                ValueError("Base graph must have at least two edges")
        else:
            print("Base structure to gmm must be a NetworkX Graph or DiGraph object, no change made.")
            
    def revert_base(self):
        """Reverts base graph to initial structure"""
        self.base=copy.deepcopy(self.original)
    
    def set_termination(self, T):
        """Set the termination rule"""
        try:
            T(self.base)
            self.termination=T
        except TypeError:
            raise TypeError("T must be a NetworkX compatible function, no change made.")
            
    def apply_termination(self):
        """Applies the termination rule to the current base graph"""
        return self.termination(self.base)
    
    def set_rule(self, R):
        """Set growth function"""
        try:
            R(self.base,self.test_graph)
            self.rule=R
        except TypeError:
            print("R must be a function compatible with NetworkX graph objects, no change made.")
            
    def apply_rule(self,new,set_result=False):
        """Applies the growth rule to the current base graph with some new structure. If set_result is
        True then set base graph to result of rule application"""
        # Graph types must match, do coercion step if necessary
        try:
            if self.base.is_directed() is True and new.is_directed() is False:
                new.to_directed()
            else:
                if new.is_directed() is True:
                    new.to_undirected()
        except TypeError:
            raise TypeError("New graph structure not a NetworkX Graph or DiGraph object")
        # Apply rule
        if set_result is True:
            self.base=self.rule(self.base,new)
            return self.base
        else:
            return self.rule(self.base,new)
            
    def am_gmm(self):
        """Simple function to test if object is a gmm"""
        return self.am_gmm
            

if __name__ == '__main__':
    # Create most basic GMM object with five node cycle graph as base.
    
    import gmm
    
    G=nx.cycle_graph(5)
    model=gmm.gmm(G)
    
    # Using five node cycle graph, create gmm object with node ceiling
    # termination rule of 100 nodes
    
    def degree_ceiling(G):
        if G.number_of_nodes()>=100:
            return True
        else:
            return False
    
    model=gmm.gmm(G,degree_ceiling)
    
    # Finally, add a simple random growth rule
    
    def rand_add(base, new):
        from numpy.random import randint
        new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
        new_base=nx.compose(base,new)
        base_connector=randint(base.number_of_nodes())
        new_connector=randint(min(new.nodes()),max(new.nodes()))
        new_base.add_edge(base_connector,new_connector)
        return new_base
    
    model.set_rule(rand_add)
