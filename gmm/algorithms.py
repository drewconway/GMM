#!/usr/bin/env python
# encoding: utf-8
"""
algorithms.py

Purpose:  Algorithms used to simulate graph structues using Graph Motif Models

Author:   Drew Conway
Email:    drew.conway@nyu.edu

"""
__author__="Drew Conway (drew.conway@nyu.edu)"
__all__=["simulate","draw_structure","motif_counts","get_motifs","all_graphs","poisson_mass"]
__docformat__ = "restructuredtext en"

import copy
import networkx as nx
from scipy import random,stats
from numpy import mean

def simulate(gmm,tau,poisson=True,new_name="GMM Simulation"):
    """
    The primary function for generating networks using the graph motif modeling technique.  
    The function takes two arguments, a gmm object and a tau value, and returns a NetworkX 
    graph object based on the construct of the gmm passed.
    
    Parameters
    ----------
    gmm : A graph motif model object.
    
    tau : An number greater than or equal to two, which is used to generate a set of graph 
        motifs for the simulation.
    
    Returns
    ----------
    gmm_sim : A NetworkX graph object derived from a graph motif model simualtion on a 
        given gmm object.
    """
    # First check that the user arguments pass inspection
    try:
        gmm.am_gmm
    except ValueError:
        raise Value("Argument to function must be well-defined gmm object")
    else:
        if tau<2:
            raise ValueError("The value for tau must be greater than or equal to two")
        else:
            # Do simulation
            while gmm.apply_termination():
                motif_dist=motif_counts(gmm,tau)    # Raw motif counts from gmm base graph
                # Poission PMF used to estimate mass for all motifs? (default)
                if poisson:
                    motif_mass=poisson_mass(motif_dist)
                # Otherwise, use count ratios
                else:
                    total_counts=sum([(c) for (a,b,c) in motif_dist])
                    motif_mass=[(a,b,float(c)/total_counts) for (a,b,c) in motif_dist]
                new_structure=draw_structure(motif_mass)
                gmm.apply_rule(new_structure,set_result=True)
        # Reset name
        gmm.get_base().name=new_name
                
                
def draw_structure(motif_mass):
    """
    Take a list of tuples of the construct (index,motif,probability mass) and takes a random
    draw of a motif based on the probability masses.
    
    Parameters
    ----------
    motf_mass : A list of tuples of the constructure (index,motif,probability mass), likely constructred
        by the gmm_simulate function.
                
    Returns 
    ----------
    motif : A randomly drawn graph motif, as a NetworkX graph object
    """
    # Truncate probabilities around tau motifs
    mass_sum=sum([(c) for (a,b,c) in motif_mass])
    probabilities=map(lambda x: motif_mass[x][2]/mass_sum,range(len(motif_mass)))

    # Step through probability mass to find appropriate 
    # motif given a uniform draw
    draw=random.uniform()
    motif_index=0
    mass_sum=probabilities[motif_index]
    while draw>mass_sum:
        motif_index+=1
        mass_sum+=probabilities[motif_index]
    return motif_mass[motif_index][1] # Return the appropriate motif


def motif_counts(gmm,tau):
    """
    Returns dictionary keyed by graph motifs and values as the number of subgraph isomorphisms 
    for the given motif counted in the base structure of the given GMM object.

    Parameters
    ----------
    gmm : A graph motif model object.

    tau : An integer greater than or equal to 2, which designates the number of nodes in the 
        largest graph in set of graph motifs used in the given model.
        
    Returns
    ----------
    subgraph_counts : A list of tuples with the following construction (index,motif,count), where 
        "count" is the number of subgraph isomorphisms for the given motif counted in the base 
        structure of the given GMM object and index is the motif index.
    """
    base=gmm.get_base()
    base_direction=base.is_directed()   # Check if GMM base is directed, motifs must match
    motif_counts=get_motifs(tau,base_direction)
    # Performing the counting of subgraph isomorphism for every motif given the base structure
    for motif in motif_counts:
        index=motif[0]
        # Calculate subgraph ismorphisms
        if base_direction:
            GM=nx.DiGraphMatcher(base,motif[1])
        else:
            GM=nx.GraphMatcher(base,motif[1])
        # Count subgraph isomorphisms and update
        count=0
        for i in GM.subgraph_isomorphisms_iter():
            count+=1
        motif_counts[index]=(index,motif_counts[index][1],count)
    return motif_counts


def get_motifs(tau,directed_motifs):
    """
    Returns a list of tupples of all possible single component non-singleton subgraphs from a dyad to 
    $G=\{V=tau,E=\frac{tau(tau-1)}{2}\}$, i.e. the K_{tau} complete graph.

    The list represents the set of graph moitifs used in a given model.  This function is 
    primarily meant as a helper to motif_counts, but can be used for other user purposes.

    Parameters
    ----------
    tau : An integer greater than or equal to 2, which designates the number of nodes in the 
        largest graph in set of graph motifs used in the given model.
        
    directed_motifs : A boolean designating whether the motifs should be directed.
        
    Returns
    ----------
    motifs : A list of tuples of the following construction (index,motif,0).  The final zero is used as 
        is a placeholder and will set the counts of subgraph isomorphism for each motif in the GMM's 
        base structure.
    """
    motifs=list()
    motif_index=0
    for v in xrange(2,tau+1):
        graphs=all_graphs(v,directed_motifs)
        for g in graphs:
            motifs.append((motif_index,g,0))
            motif_index+=1
    return motifs
        
        
def all_graphs(num_nodes,directed_motifs):
    """
    Retuns a list of all possible single component graphs given some number of nodes. This function 
    is meant primarily as a helper function to get_subgraphs, but can be used for other user purposes.

    Parameters
    ----------
    num_nodes : The number of nodes in the complete graph

    directed_motifs : A boolean designating whether the motifs should be directed.

    Returns
    ----------
    graphs : A list of all possible single component graphs given some number of nodes
    """
    graphs=list()
    complete=nx.complete_graph(num_nodes)    # Start with complete graph
    complete_copy=copy.deepcopy(complete)
    if directed_motifs:
        complete=complete.to_directed()
        complete_copy=complete_copy.to_directed()
        complete_ud=all_graphs(num_nodes,False) # RECURSION, FTW!
    edges=complete.edges()
    while complete.number_of_edges()>0:
        # Iteratively remove edges, and capture single component subgraphs
        e=edges.pop()
        complete.remove_edge(e[0],e[1])
        if directed_motifs:
            if nx.number_weakly_connected_components(complete)==1:
                graphs.append(copy.deepcopy(complete))
        else:
            if nx.number_connected_components(complete)==1:
                graphs.append(copy.deepcopy(complete))
    graphs.append(complete_copy)
    # Add recursively produced undirected graphs to directed set
    if directed_motifs:
        graphs_edges=map(lambda g: g.edges(), graphs)
        graphs_edges.sort()
        for i in complete_ud:
            ud_edges=i.edges()
            ud_edges.sort()
            if ud_edges not in graphs_edges:
                graphs.append(i.to_directed())
    return graphs
    
    
def poisson_mass(motif_counts):
    """
    Will calculate a shape parameter (lambda) for a Possion distribution given some motif counts.
    
    Parameters
    ----------
    motif_counts :  A dictionary of motif counts keyed by (order, motif) tuples with count values.
    
    Returns
    ----------
    motif_pmf : A dictionary of motif counts keyed by (order, motif) tuples with Poisson probability
        mass function as values given the lambda estimate from the motif_counts.
    
    Notes
    -----
    There are two critical assumption inherent to using the Poisson fit method for estimating:
    
    - Graphs can be ordered by complexity.  This is tacitly implied in the model, but when using
    the Poisson fit to estimate mass for all motifs it is explicit.
    
    - Perhaps more importantly, motif counts are distributed Possion
        
    In this case the most basic motif, the dyad, maps to k=0, and the assumed k value increases as 
    both node and edge complexity increase in motifs.
        
    The primary advantage of this method is to assure that there is some positive probability mass
    for every motif in the model. In some cases a motif may have a zero probability given some
    base structure.  This method will prevent that outcome.
    """
    counts=[(c) for (a,b,c) in motif_counts]
    estimate=mean(counts)
    # Use SciPy's Poisson PMF calculation function
    motif_pmf=map(lambda x: stats.poisson.pmf(x,estimate),range(len(counts)))
    for i in range(len(motif_pmf)):
        motif_counts[i]=(motif_counts[i][0],motif_counts[i][1],motif_pmf[i])
    return(motif_counts)

        

if __name__ == '__main__':
    pass
    