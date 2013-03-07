#!/usr/bin/env python
# encoding: utf-8
"""
actor_projection.py

Purpose:  Create bipartite projections that retain attribute data from original graph

Author:   Drew Conway
Email:    drew.conway@nyu.edu
Date:     2010-11-11

Copyright (c) 2010, under the Simplified BSD License.  
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.
"""

import sys
import os
import networkx as nx


def main():
    # Load data an symmetrize
    ssrn_mc=nx.read_graphml("../data/ssrn_mc.graphml")
    ssrn_mc=ssrn_mc.to_undirected()
    
    # Create projections
    author_graph=nx.project(ssrn_mc,[(a) for (a,b) in ssrn_mc.nodes(data=True) if b["type"]=="1"])
    article_graph=nx.project(ssrn_mc,[(a) for (a,b) in ssrn_mc.nodes(data=True) if b["type"]=="0"])
    
    # Reset vertex attibutes
    author_graph.add_nodes_from([(a,b) for (a,b) in ssrn_mc.nodes(data=True) if b["type"]=="1"])
    article_graph.add_nodes_from([(a,b) for (a,b) in ssrn_mc.nodes(data=True) if b["type"]=="0"])
    
    # Write graph data
    nx.write_graphml(author_graph,"../data/authors_mc.graphml")
    nx.write_graphml(article_graph,"../data/articles_mc.graphml")
    
if __name__ == '__main__':
	main()

