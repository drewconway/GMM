# File-Name:       cluster_analysis.R           
# Date:            2010-11-12                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Perform basic community detection on SSRN affilaitions graphs and
#                   add partition data to each vertex, fix labels, and re-save.
# Data Used:       authors_mc.graphml,articles_mc.graphml,abstracts.csv
# Packages Used:   igraph
# Output File:    
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2010, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

# Load libraries and graph data, and remove id attribute
library(igraph)
authors.graph<-read.graph("../data/authors_mc.graphml", format="graphml")
articles.graph<-read.graph("../data/articles_mc.graphml", format="graphml")
authors.graph<-set.vertex.attribute(authors.graph, name="label", value=V(authors.graph)$name)
articles.graph<-set.vertex.attribute(articles.graph, name="label", value=V(articles.graph)$name)
authors.graph<-remove.vertex.attribute(authors.graph,"id")
authors.graph<-remove.vertex.attribute(authors.graph,"name")
articles.graph<-remove.vertex.attribute(articles.graph,"id")
articles.graph<-remove.vertex.attribute(articles.graph,"name")

# Perform three common community detection alogorithms available in igrapgh

# Add spinglass community
authors.spinglass<-spinglass.community(authors.graph)
articles.spinglass<-spinglass.community(articles.graph)
authors.graph<-set.vertex.attribute(authors.graph, name="spinglass", value=as.character(authors.spinglass$membership))
articles.graph<-set.vertex.attribute(articles.graph, name="spinglass", value=as.character(articles.spinglass$membership))

# Add walktrap community
authors.walktrap<-walktrap.community(authors.graph, weight=NULL)
articles.walktrap<-walktrap.community(articles.graph, weight=NULL)
authors.graph<-set.vertex.attribute(authors.graph, name="walktrap", value=as.chatacter(authors.walktrap$membership))
articles.graph<-set.vertex.attribute(articles.graph, name="walktrap", value=as.character(articles.walktrap$membership))

# Add fastgreedy community
authors.fastgreedy<-fastgreedy.community(authors.graph)
articles.fastgreedy<-fastgreedy.community(articles.graph)
authors.fastgreedy.membership<-community.to.membership(authors.graph,merges=authors.fastgreedy$merges,which.max(authors.fastgreedy$modularity))
articles.fastgreedy.membership<-community.to.membership(articles.graph,merges=articles.fastgreedy$merges,which.max(articles.fastgreedy$modularity))
authors.graph<-set.vertex.attribute(authors.graph, name="fastgreedy", value=as.character(authors.fastgreedy.membership$membership))
articles.graph<-set.vertex.attribute(articles.graph, name="fastgreedy", value=as.character(articles.fastgreedy.membership$membership))

# Save new graphs 
write.graph(authors.graph, "../data/authors_community.graphml", format="graphml")
write.graph(articles.graph, "../data/articles_community.graphml", format="graphml")

# Add article community information to abstract data for topic modeling

abstracts<-read.csv("../data/abstracts.csv", stringsAsFactors=FALSE)
abstracts<-abstracts[match(V(articles.graph)$label,abstracts$Article),]
communities<-as.data.frame(cbind(V(articles.graph)$label,articles.spinglass$membership,articles.walktrap$membership,articles.fastgreedy.membership$membership),stringsAsFactors=FALSE)
names(communities)<-c("Article","Spinglass","Walktrap","Fastgreedy")

# Sort for merging
abstracts<-abstracts[with(abstracts, order(Article)),]
communities<-communities[with(communities, order(Article)), ]
abstracts.community<-merge(abstracts, communities, by.all=Article)

# Re-write abstract data
write.csv(abstracts.community, "../data/abstracts_communities.csv", row.names=FALSE)