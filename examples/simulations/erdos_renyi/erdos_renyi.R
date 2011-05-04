# File-Name:       erdos_renyi.R           
# Date:            2011-04-26                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Calculate goodness of fit between degree distributions
#                   for ER and GMM models
# Data Used:       data/*
# Packages Used:   igraph 
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2011, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.   

library(igraph)
library(vcd)
library(ggplot2)                                                     

degree.dist<-read.csv("data/degree_dist.csv", stringsAsFactors=TRUE)
data.dir<-"data/"
graph.files<-dir(data.dir)[-1]

# A function to calculate a GLM with binomial link
binom.fit<-function(path) {
    G<-read.graph(paste(data.dir,path,sep=""), format="graphml")
    dd<-degree.distribution(G)
    binom.dens<-dbinom(0:max(degree(G)), length(V(G))-1, 0.5)
    return(cbind(path,0:max(degree(G)),dd,binom.dens))
}

# Collect all of the fit data into a single data frame
graph.fits<-lapply(graph.files, binom.fit)
fit.df<-data.frame(do.call(rbind, graph.fits), stringsAsFactors=FALSE)
names(fit.df)<-c("graph","degree","deg.dist","binom")
fit.df$graph<-as.factor(fit.df$graph)
fit.df$degree<-as.numeric(fit.df$degree)
fit.df$deg.dist<-as.numeric(fit.df$deg.dist)
fit.df$binom<-as.numeric(fit.df$binom)

for(g in levels(fit.df$graph)) {
    graph.df<-subset(fit.df, graph==g)
    fit.plot<-ggplot(graph.df, aes(xmin=degree-.5, xmax=degree+.5, ymin=0, ymax=deg.dist))+
        geom_rect(aes(fill="darkgrey"))+
        geom_point(aes(x=degree, y=binom, color="darkred"))+
        geom_line(aes(x=degree, y=binom, color="darkred"))+
        scale_fill_manual(values=c("darkgrey"="darkgrey"), name="Observed Degree", breaks="darkgrey", label="")+
        scale_color_manual(values=c("darkred"="darkred"), name="Fitted Binomial", breaks="darkred", label="")+
        scale_x_continuous(limits=c(0,65))+xlab("Degree")+ylab("Density")+theme_bw()
    ggsave(plot=fit.plot, filename=paste("images/",as.character(g),"_fitplot.pdf",sep=""))
}

# Fit distributions using minimum Chi-squared
fit.lm<-function(path) {
    G<-read.graph(paste(data.dir,path,sep=""), format="graphml")
    deg<-degree.distribution(G)
    g.bn<-dbinom(0:max(degree(G)), size=length(V(G))-1, prob=0.5)
    return(lm(deg ~ g.bn))
    
}   

graph.names<-levels(fit.df$graph)
graph.lm<-lapply(graph.files, fit.lm)
graph.rmse<-sapply(graph.lm, function(lm) sqrt(mean(lm$residual^2)))
graph.aic<-sapply(graph.lm, AIC)
names(graph.chisqr)<-graph.names



