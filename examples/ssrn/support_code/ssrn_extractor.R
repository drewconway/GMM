# File-Name:       ssrn_extractor.R           
# Date:            2010-11-09                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         To build a citation network from the SSRN Conflict Studies eJournal
# Data Used:       
# Packages Used:   XML, igraph, plyr
# Output File:    
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2010, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

require(XML)
require(igraph)
library(plyr)

#### STEP 1: Get initial base abstract IDs for all artciles in SSRN Conflict Studies eJournal

# Must hard-code article search numbers at this point
# This max artcile number is only valid as of 2010-12-11, go to the link below update for your needs
# http://papers.ssrn.com/sol3/JELJOUR_Results.cfm?form_name=journalBrowse&journal_id=998694
max.articles<-2554
article.seq<-seq(1,max.articles,20)

# SSRN search URL constants, bound by search sequence
url.front<-"http://papers.ssrn.com/sol3/Jeljour_results.cfm?nxtres="
url.back<-"&form_name=journalBrowse&journal_id=998694&Network=no&SortOrder=ab_approval_date&stype=desc&lim=false"

# Define function to extract unique abstract ID each artcile on a given search result page
get.abs<-function(ssrn.url) {
    # SSRN's html is grossly malformed. Lots of text manipulation magic happening here.
    doc<-htmlTreeParse(ssrn.url)
    raw.html<-as.vector(unlist(doc$children$html[["body"]][["table"]][[3]][[3]]))
    abs<-raw.html[which(grepl("abstract_id",raw.html,fixed=TRUE))]
    return(do.call("rbind",strsplit(abs,"=",fixed=TRUE))[,2])
}

# Now get all abstracts
# WARNING: this can take several minutes, only run if abstract ID have not yet been parsed
all.abs<-unlist(lapply(article.seq, function(x) { get.abs(paste(url.front,x,url.back,sep="")) }))

#### STEP 2: Create the network from above seeds

# SSRN abstract URL constants
abs.url<-"http://papers.ssrn.com/sol3/papers.cfm?abstract_id="

# Define a function that takes an SSRN absract URL and extracts author IDs
get.doc<-function(abs.url) {
    # Same issues with ugly HTML
    doc<-htmlTreeParse(abs.url)
    raw.html<-as.vector(unlist(doc$children$html[["body"]][[6]][[1]]))
    # Get authors
    authors<-raw.html[which(grepl("per_id",raw.html,fixed=TRUE))]
    authors<-do.call(rbind, strsplit(authors,"=",fixed=TRUE))[,2]
    # Get abstract, and use as edge data
    abs.local<-which(grepl("Abstract:",raw.html,fixed=TRUE))
    abs<-raw.html[abs.local+8]
    date.local<-which(grepl("Date posted:",raw.html,fixed=TRUE))
    date<-raw.html[date.local+4]
    if(length(authors)<1) {
        return(list(Authors=cbind(NA,date),Abstract=abs))
    }
    else {
        return(list(Authors=cbind(authors,date), Abstract=abs))
    }
}
    

# Define function that takes a abstract ID, and returns co-author edgelist
doc.el<-function(abs.id) {
    doc.info<-get.doc(paste(abs.url,abs.id,sep=""))
    return(list(EL=cbind(doc.info[["Authors"]],abs.id), Abstract=doc.info[["Abstract"]]))
}

# Build full network. WARNING: this takes A LOT longer than the initial abstract parsing.
# Sometimes over an hour, so grab some coffee and relax.
ssrn.data<-lapply(all.abs,doc.el)

# Create a single two column matrix, then convert to data frame
ssrn.el<-lapply(ssrn.data, function(d) d[["EL"]])
ssrn.df<-do.call("rbind",ssrn.el)

# Shift Article to second column and convert dates to POSIX
ssrn.df<-as.data.frame(list(Author=ssrn.df[,1],Article=ssrn.df[,3],Date=as.Date(ssrn.df[,2], "%B %d, %Y")), stringsAsFactors=FALSE)

# Clean out duplicates, but retain dup counts as edge weight,
# though not really a weight.
ssrn.clean<-ddply(ssrn.df,.(Author,Article), nrow)
names(ssrn.clean)<-c("Author","Article","weight")

# Create data frame of just Artcile IDs and Abstract for future topic modeling
abs.id<-function(ssrn.list) {
    # Function returns [article.id,abs] vector
    article.id<-ssrn.list[["EL"]][1,2]
    abs<-ssrn.list[["Abstract"]]
    return(cbind(article.id,abs))
}
abs.df<-lapply(ssrn.data, abs.id)
abs.df<-do.call("rbind", abs.df)
abs.df<-as.data.frame(abs.df, stringsAsFactors=FALSE)
names(abs.df)<-c("Article", "Abstract")

# Save data frames
write.csv(ssrn.df, "../data/raw_data.csv", row.names=FALSE)
write.csv(ssrn.clean, "../data/clean_data.csv", row.names=FALSE)
write.csv(abs.df, "../data/abstracts.csv", row.names=FALSE)

#### STEP 3: Create networks and save

ssrn.df<-read.csv("../data/raw_data.csv", stringsAsFactors=FALSE)
ssrn.clean<-read.csv("../data/clean_data.csv", stringsAsFactors=FALSE)

# Create igraph objects, and type attribute, and save as GraphML

# Temporal network
temporal.graph<-graph.data.frame(ssrn.df)
temporal.type<-is.bipartite(temporal.graph)$type
temporal.type<-ifelse(temporal.type,0,1)
temporal.graph<-set.vertex.attribute(temporal.graph, name="type", value=as.character(temporal.type))
write.graph(temporal.graph,"../data/ssrn_temporal.graphml",format="graphml")


# Weighted by counts
weighted.graph<-graph.data.frame(ssrn.clean)
weighted.type<-is.bipartite(weighted.graph)$type
weighted.type<-ifelse(weighted.type,0,1)
weighted.graph<-set.vertex.attribute(weighted.graph, name="type", value=as.character(weighted.type))
write.graph(weighted.graph,"../data/ssrn_weighted.graphml",format="graphml")

# Extract main component and save
graph.components<-decompose.graph(weighted.graph,mode="weak", min.vertices=3)   # ignore single author dyads
component.sizes<-sapply(1:length(graph.components), function(x) vcount(graph.components[[x]]))
main.component<-graph.components[[which.max(component.sizes)]]
write.graph(main.component, "../data/ssrn_mc.graphml",format="graphml")