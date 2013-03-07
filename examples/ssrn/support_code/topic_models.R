# File-Name:       topic_models.R           
# Date:            2010-11-13                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Generate topic models from artcile abstracts by community with LDA 
# Data Used:       abstracts_community.csv
# Packages Used:   topicmodels
# Output File:     topic_models.csv
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2010, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

# Load libraries and data
library(topicmodels)
library(reshape)
library(plyr)
abstracts<-read.csv("../data/abstracts_communities.csv", stringsAsFactors=FALSE)

# Create corpus from Abstracts for all community factors
create.lda<-function(abstracts, communities, num.topics=10, num.words=10, stemming=FALSE, stopwords = TRUE, minWordLength = 4) {
    abs.matrix<-cbind(abstracts,communities)
    # Create list for LDA terms
    ssrn.terms<-list()
    com.levels<-unique(communities)
    dtm.control<-list(stemming = stemming, stopwords = stopwords, minWordLength = minWordLength)
    for(c in com.levels) {
        com.corp<-Corpus(VectorSource(abs.matrix[which(communities==c),1]))
        dtm<-DocumentTermMatrix(com.corp,control=dtm.control)
        dtm<-removeSparseTerms(dtm,0.95)
        ssrn.terms[[as.character(c)]]<-tryCatch(get_terms(LDA(dtm,control = list(alpha = 0.05), k = num.words),num.topics), error=function(e) NA)
    }
    return(ssrn.terms)
}

# For some reason community 4 in the Fastgreedy partition is causing an error in
# the LDA function, so at this point that group is set to null.
topic.models<-lapply(3:5, function(c) create.lda(abstracts$Abstract, abstracts[,c], num.topics=5, num.words=5))
names(topic.models)<-c("Spinglass","Walktrap","Fastgreedy")

# Create data frame from list
topic.df<-melt(topic.models)
names(topic.df)<-c("Term.Num","Topic", "Term", "Partition", "Community")
topic.df$Term<-as.character(topic.df$Term)
topic.df$Term<-gsub("[[:punct:]]", "", topic.df$Term)  # Strip out punctuation

# Create intra-model/partion counts
topic.counts<-ddply(topic.df,.(Term,Partition,Community), nrow)
names(topic.counts)<-c("Term","Partition","Community", "Count")

# Output topic models
write.csv(topic.counts, "../data/topic_models.csv", row.names=FALSE)
