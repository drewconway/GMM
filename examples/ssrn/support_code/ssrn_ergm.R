# File-Name:       ssrn_ergm.R           
# Date:            2011-01-19                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Fit and simulate ERGM model from SSRN 2009 network
# Data Used:       ssrn_2009.graphml
# Packages Used:   statnet
# Output File:    
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2011, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

# Load libraries and data
library(statnet)

ssrn.data<-as.matrix(read.table("../data/ssrn_2009.txt"))
ssrn.net<-network(ssrn.data,bipartite=TRUE)

# Specify various ERGM models for the netowrk
ssrn.ergm1<-ergm(ssrn.net ~ degree(1:2), burnin=80000)
ssrn.ergm2<-ergm(ssrn.net ~ degree(1:2) + kstar(1:2), burnin=80000)
ssrn.ergm3<-ergm(ssrn.net ~ degree(1:2) + kstar(1:2) + density, burnin=80000)
ssrn.ergm4<-ergm(ssrn.net ~ degree(1:2) + kstar(1:2) + density + gwdegree(1), burnin=80000)

# Simulate networks from models
ssrn.sim1<-simulate(ssrn.ergm1)
ssrn.sim2<-simulate(ssrn.ergm2)
ssrn.sim3<-simulate(ssrn.ergm3)
ssrn.sim4<-simulate(ssrn.ergm4)

sim1.el<-as.edgelist(ssrn.sim1)
sim2.el<-as.edgelist(ssrn.sim2)
sim3.el<-as.edgelist(ssrn.sim3)
sim4.el<-as.edgelist(ssrn.sim4)

# Load igraph, and save simulation as GraphML
library(igraph)

write.graph(graph(sim1.el),"../data/ssrn_sim1.graphml",format="graphml")
write.graph(graph(sim2.el),"../data/ssrn_sim2.graphml",format="graphml")
write.graph(graph(sim3.el),"../data/ssrn_sim3.graphml",format="graphml")
write.graph(graph(sim4.el),"../data/ssrn_sim4.graphml",format="graphml")