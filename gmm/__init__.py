#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Purpose:  High-level load for graph motif model classes and functions

Author:   Drew Conway
Email:    drew.conway@nyu.edu

"""

import sys
import os

# Required third-party 
try:
    import networkx as nx
except ImportError:
    raise ImportError("You must install NetworkX (http://networkx.lanl.gov/) to use gmm")
try:
    from scipy import random,stats
    from numpy import mean
except ImportError:
    print("You must install SciPy (http://www.scipy.org/) to use gmm")
    
# GMM classes and functions
import gmm
from gmm import *
import algorithms
from algorithms import *
