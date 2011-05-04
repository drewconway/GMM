#!/usr/bin/env python
# encoding: utf-8
"""
barabasi_albert.py

Purpose:  Recover the Barabasi-Albert "preferrential attachment" random graph model with GMM
          see, A. L. Barabási and R. Albert "Emergence of scaling in random networks", Science 286, pp 509-512, 1999.

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2011-03-28

Copyright (c) 2011, under the Simplified BSD License.  
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""
import sys
import os
import networkx as nx
import gmm
from numpy import arange, random, mean
import matplotlib.pylab as plt
import smtplib
from datetime import datetime
        

def barabasi_albert_simulation(test_graph, m, seed=None, verbose=False):
    """
    Given a set of Barabasi-Albert "preferential attachment" random graphs, simulate using above growth rule
    """
    
    # Required to change the growth rule dynamically, so 
    # we define it inside the simulation function
    def barabasi_albert_growth(base, new):
        """
        Select m random nodes from new_nodes and connect each node to the nodes in base_nodes
        as a function of the degree of the node in base_nodes.  The basic "preferential 
        attachment" model.
        """

        # To keep new nodes from over-writing current ones rename the new nodes starting
        # from the last node in base
        new=nx.convert_node_labels_to_integers(new,first_label=max(base.nodes())+1)
        new_nodes=new.nodes()
        base_nodes=base.nodes()
        new_base=nx.compose(base,new)

        # Shuffle new_nodes
        random.shuffle(new_nodes)

        # Create edge test based on degree centrality
        base_degree=nx.degree_centrality(base).items()

        # Create "preferential attachment" structure by connecting m
        # nodes from new structure to nodes in base as a function of 
        # base nodes' degree centrality
        for i in xrange(m):
            edge_made=False
            while edge_made is False:
                # Randomly select a node in base and add connection
                # based on its degree centrality
                p=random.uniform()
                j=random.randint(0,len(base_nodes))
                if p <= base_degree[j][1]:
                    k=random.randint(len(new_nodes))    # Randomly select a new node
                    new_base.add_edge(new_nodes[k],base_degree[j][0])
                    edge_made=True
        return(new_base)
    
    # Simple node ceiling of 1,000 nodes 
    def node_ceiling(G):
        if G.number_of_nodes()>=1000:
            return False
        else:
            return True    
        
    # Setup GMM
    barabasi_albert_model=gmm.gmm(test_graph)
    barabasi_albert_model.set_termination(node_ceiling)
    barabasi_albert_model.set_rule(barabasi_albert_growth)
    sim_name="GMM SIMULATION-"+barabasi_albert_model.get_base().name
    gmm.algorithms.simulate(barabasi_albert_model,tau=3, poisson=True, seed=seed, new_name=sim_name)
    
    # Run simulation
    simualted_barabasi_albert=barabasi_albert_model.get_base()
    
    # Print the results to stdout
    if verbose:
        print(nx.info(simualted_barabasi_albert))
        print("")
        
    # Return a list of simulated graphs
    return(simualted_barabasi_albert)
    
def noticeEMail(starttime):
    """Sends an email message once the script is completed"""
    runtime=datetime.now() - starttime
    
    fromaddr='agconway@gmail.com'
    toaddr='agconway@gmail.com'
    
    server=smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login('agconway','diami2435343')
    
    senddate=datetime.strftime(datetime.now(), '%Y-%m-%d')
    subject="Your AWS job has completed"
    m="Date: %s\r\nFrom: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\n\r\n" % (senddate, fromaddr, toaddr, subject)
    msg='''
    
    Job runtime: '''+str(runtime)
    
    server.sendmail(fromaddr, toaddr, m+msg)

def main():
    # Start time of script for email
    starttime=datetime.now()
    
    """
    1. Generate a set of Barabasi-Albert graphs to experiment on
    """

    # Set seed for random graph
    random_seed=851982
    random.seed(random_seed)    # Seed the generator
    
    # Test ranges from original Barabasi-Albert paper
    m_range=[1,3,5,7]
    n_ba=1000

    # Produce a set of Barabasi-Albert random graphs with 300,000 nodes and m values ranging
    # from 1, 3, 5, 7.  These values were chosen as they match those used in the original
    # BA paper.
    test_set=map(lambda i: nx.barabasi_albert_graph(n_ba, m=i), m_range)

    # Output the base graphs as GraphML
    for g in test_set:
       nx.write_graphml(g, "data/base/"+g.name+"_"+str(g.number_of_nodes())+".graphml") 

    print("Base Barabasi-Albert graphs generated")
    
    """
    2. We now perform a grid-search over the space for the number of connections to make in the
    Barabasi-Alebrt model (m) and the size of the base structure from which the GMM will simlate a 
    Barabasi-Alebrt netwotk using the above growth rule.
    
    Barabasi-Alebrt uses the scaling parameter of the degree distributions to test their model, so
    this metrics will be used to compare the GMM simualtions versus the networks produced using the 
    Barabasi-Alebrt model.
    """
    
    # Variation in sizes of base
    base_sizes=arange(100,1000,100)
    
    # Run simualtions and output
    graph_id=0
    for x in range(28): # Runs this a bunch of time to produce a 
        for m in m_range:
            for n in base_sizes:
                test_graph=nx.barabasi_albert_graph(n, m)
                print("id = "+str(graph_id)+"\nm = "+str(m)+"\nn = "+str(n))+"\n\n"
                gmm_simulation=barabasi_albert_simulation(test_graph, m, seed=random_seed, verbose=False)
                nx.write_graphml(gmm_simulation, "data/simulated/"+str(graph_id)+"_barabas_albert_graph_"+str(m)+"_"+str(n)+"_.graphml")
                graph_id+=1
                del gmm_simulation
                del test_graph
                
    
    # Send email to notify me we're done
    noticeEMail(starttime)


if __name__ == '__main__':
	main()

