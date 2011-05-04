# File-Name:       ba_analysis.R           
# Date:            2011-04-07                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Statistically explore the results of the Barabasi-Albert GMM models
# Data Used:       Output from barabasi_albert.py
# Packages Used:   igraph, ggplot2
# Output File:    
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2011, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

# Load libraries
library(igraph)
library(ggplot2)
library(stats4)

# A function that retuns the scaling parameter for the degree 
# distribution of a given igraph network
get.scaling<-function(graph) {
    d<-degree(graph)
    x<-0:max(d)
    dd<-c(0,hist(d, plot=FALSE, breaks=x)$counts)
    dd.log<-ifelse(dd<1,0,log10(dd))
    x.log<-ifelse(x<1,0,log10(x))
    dd.fit<-lm(dd.log ~ x.log)
    return(as.numeric(dd.fit$coeff)[2])
}

# A function that returns a vector of c(m, base_size, scaling) for
# some graphml file produced by a Barabasi-Albert GMM simulation
ba.test<-function(file) {
    g<-read.graph(file, format="graphml")
    params<-strsplit(file,"[_.]")
    return(c(params[[1]][5],params[[1]][6],get.scaling(g)))
}

base.test<-function(file) {
    g<-read.graph(file, format="graphml")
    params<-strsplit(file,"[_.]")
    return(c(params[[1]][2],get.scaling(g)))
}

get.plfit<-function(file, base=FALSE) {
    g<-read.graph(file, format="graphml")
    deg<-degree(g)
    graph.name<-strsplit(file, "[/_.]")
    if(base) {
        return(c(graph.name[[1]][3],graph.name[[1]][4],coef(power.law.fit(deg))))
    }
    else {
        return(c(graph.name[[1]][7],graph.name[[1]][8],coef(power.law.fit(deg))))
    }
}


# All network files
sim.dir<-"data/simulated/"
sim.files<-names(file.access(dir(sim.dir)))

base.dir<-"data/base/"
base.files<-names(file.access(dir(base.dir)))

# Get scaling data for all simulations
scaling.data<-lapply(sim.files, function(f) ba.test(paste(sim.dir,f,sep="")))
scaling.matrix<-do.call(rbind, scaling.data)
scaling.df<-data.frame(scaling.matrix, stringsAsFactors=FALSE)

# Clean up
names(scaling.df)<-c("m","base.size","scaling")
scaling.df$m<-as.factor(scaling.df$m)
scaling.df$base.size<-as.numeric(scaling.df$base.size)
scaling.df$scaling<-abs(as.numeric(scaling.df$scaling))
scaling.df<-scaling.df[with(scaling.df, order(m,base.size)),]


# Now do the same for a bunch of simulated BA networks
base.data<-lapply(base.files, function(f) base.test(paste(base.dir,f,sep="")))
base.matrix<-do.call(rbind, base.data)
base.df<-data.frame(base.matrix, stringsAsFactors=FALSE)

# Clean up
names(base.df)<-c("m","scaling")
base.df$m<-as.factor(base.df$m)
base.df$scaling<-abs(as.numeric(base.df$scaling))
base.df<-base.df[with(base.df, order(m)),]


# Next, use MLE to estimate the PL fit
alpha.sims<-lapply(sim.files, function(f) get.plfit(paste(sim.dir,f,sep="")))
alpha.sims.matrix<-do.call(rbind, alpha.sims)
alpha.sims.df<-data.frame(alpha.sims.matrix, stringsAsFactors=FALSE)
names(alpha.sims.df)<-c("m", "base.size", "alpha")
alpha.sims.df$m<-as.factor(alpha.sims.df$m)
alpha.sims.df$base.size<-as.numeric(alpha.sims.df$base.size)
alpha.sims.df$alpha<-as.numeric(alpha.sims.df$alpha)

alpha.base<-lapply(base.files, function(f) get.plfit(paste(base.dir,f,sep=""), base=TRUE))
alpha.base.matrix<-do.call(rbind, alpha.base)
alpha.base.df<-data.frame(alpha.base.matrix, stringsAsFactors=FALSE)
names(alpha.base.df)<-c("id", "m", "alpha")
alpha.base.df$id<-as.numeric(alpha.base.df$id)
alpha.base.df$m<-as.factor(alpha.base.df$m)
alpha.base.df$alpha<-as.numeric(alpha.base.df$alpha)

# Produce some visualizations

# m=1 plot
m1.plot<-ggplot(subset(base.df, m==1), aes(x=scaling))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(scaling.df, m==1), aes(x=scaling, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(0,2))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=1]")+
    xlab("Scaling parameter")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m1.plot, filename="m1.pdf", height=7, width=11)
    
# m=3 plot
m3.plot<-ggplot(subset(base.df, m==3), aes(x=scaling))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(scaling.df, m==3), aes(x=scaling, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(0,2))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=3]")+
    xlab("Scaling parameter")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m3.plot, filename="m3.pdf", height=7, width=11)

# m=5 plot
m5.plot<-ggplot(subset(base.df, m==5), aes(x=scaling))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(scaling.df, m==5), aes(x=scaling, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(0,2))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=5]")+
    xlab("Scaling parameter")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m5.plot, filename="m5.pdf", height=7, width=11)

# m=7 plot
m7.plot<-ggplot(subset(base.df, m==7), aes(x=scaling))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(scaling.df, m==7), aes(x=scaling, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(0,2))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=7]")+
    xlab("Scaling parameter")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m7.plot, filename="m7.pdf", height=7, width=11)

# Repeat for MLE estimates
m1.alpha.plot<-ggplot(subset(alpha.base.df, m==1), aes(x=alpha))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(alpha.sims.df, m==1), aes(x=alpha, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(1,3))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=1]")+
    xlab("MLE Power-law Alpha Estimate")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m1.alpha.plot, filename="m1_alpha.pdf", height=7, width=11)

m3.alpha.plot<-ggplot(subset(alpha.base.df, m==3), aes(x=alpha))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(alpha.sims.df, m==3), aes(x=alpha, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(1,3))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=3]")+
    xlab("MLE Power-law Alpha Estimate")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m3.alpha.plot, filename="m3_alpha.pdf", height=7, width=11)

m5.alpha.plot<-ggplot(subset(alpha.base.df, m==5), aes(x=alpha))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(alpha.sims.df, m==5), aes(x=alpha, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(1,3))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=5]")+
    xlab("MLE Power-law Alpha Estimate")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m5.alpha.plot, filename="m5_alpha.pdf", height=7, width=11)

m7.alpha.plot<-ggplot(subset(alpha.base.df, m==7), aes(x=alpha))+stat_density(aes(fill="BA", alpha=.65))+
    stat_density(data=subset(alpha.sims.df, m==7), aes(x=alpha, fill="GMM",alpha=.65))+
    scale_fill_manual(values=c("BA"="darkred", "GMM"="darkblue"), name="Model type")+
    scale_alpha(legend=FALSE)+scale_x_continuous(limits=c(1,3))+
    opts(title="Comparison between Barabasi-Albert Random\nGraph Model and GMM Equivalent [n=100, m=7]")+
    xlab("MLE Power-law Alpha Estimate")+
    ylab("Density")+
    theme_bw()
ggsave(plot=m7.alpha.plot, filename="m7_alpha.pdf", height=7, width=11)
