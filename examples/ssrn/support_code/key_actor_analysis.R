# File-Name:       key_actor_analysis.R           
# Date:            2010-11-13                                
# Author:          Drew Conway
# Email:           drew.conway@nyu.edu                                      
# Purpose:         Do basic key actor analysis on articles and authors graphs and degree
#                   and degree distribution analyses.
# Data Used:       authors_community.graphml, articles_community.graphml
# Packages Used:   igraph, ggplot2, vcd
# Output File:    
# Data Output:     
# Machine:         Drew Conway's MacBook Pro

# Copyright (c) 2010, under the Simplified BSD License.  
# For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# All rights reserved.                                                         

# Load libraries and data
library(igraph)
library(ggplot2)
library(vcd)

authors<-read.graph("../data/authors_community.graphml",format="graphml")
articles<-read.graph("../data/articles_community.graphml",format="graphml")

#### Calculate centrality measures on two-core

# Eigenvector centrality
authors.evcent<-evcent(authors)$vector
articles.evcent<-evcent(articles)$vector

# Betweenness centrality
authors.betw<-betweenness(authors)
articles.betw<-betweenness(articles)

# Closeness centrality
authors.close<-closeness(authors)
articles.close<-closeness(articles)

# Page Rank
authors.pr<-page.rank(authors)$vector
articles.pr<-page.rank(articles)$vector

# Create data frames
authors.df<-as.data.frame(list(Authors=V(authors)$label, Eig=authors.evcent, Bet=authors.betw, Close=authors.close, PR=authors.pr), stringsAsFactors=FALSE)
articles.df<-as.data.frame(list(Articles=V(articles)$label, Eig=articles.evcent, Bet=articles.betw, Close=articles.close, PR=articles.pr), stringsAsFactors=FALSE)

# Simple linear models to measure residual
authors.lm<-lm(Bet~PR,data=authors.df)
authors.res<-as.vector(authors.lm$res)

articles.lm<-lm(Bet~PR,data=articles.df)
articles.res<-as.vector(articles.lm$res)

authors.df$Res<-authors.res
articles.df$Res<-articles.res

# Create plots. Using Closeness centrality over Eigenvector because of highly localized clustering in affiliations networks
key.authors<-ggplot(subset(authors.df,Bet>0), aes(x=PR,y=Bet))+geom_text(aes(colour="ID",size=abs(Res),label=Authors,alpha=0.7))+
    theme_bw()+xlab("Page Rank")+ylab("Betweenness Centrality")+scale_size(legend=FALSE)+
    scale_colour_manual(values=c("ID"="black"),name="Authors",legend=FALSE)+scale_alpha(legend=FALSE)+
    scale_x_continuous(limits=c(min(authors.df$PR)-0.0005,max(authors.df$PR)+0.0005))+
    opts(title="Key Actors Analysis for Author Network")
ggsave(plot=key.authors, filename="../images/key_authors.pdf",height=5,width=8)
ggsave(plot=key.authors, filename="../images/key_authors.png",height=5,width=8)

key.articles<-ggplot(subset(articles.df,Bet>0), aes(x=PR,y=Bet))+geom_text(aes(colour="ID",size=abs(Res),label=Articles,alpha=0.7))+
    theme_bw()+xlab("Page Rank")+ylab("Betweenness Centrality")+scale_size(legend=FALSE)+
    scale_colour_manual(values=c("ID"="black"),name="Articles",legend=FALSE)+scale_alpha(legend=FALSE)+
    scale_x_continuous(limits=c(min(authors.df$PR)-0.0005,max(authors.df$PR)+0.0005))+
    opts(title="Key Actors Analysis for Articles Network")
ggsave(plot=key.articles, filename="../images/key_articles.pdf",height=5,width=8)
ggsave(plot=key.articles, filename="../images/key_articles.png",height=5,width=8)

#### Dergee distribution analysis

# Author network degree distribution plot
authors.degree<-degree(authors)
authors.dd<-degree.distribution(authors)
authors.table<-table(authors.degree)
authors.match<-match(0:max(authors.degree),row.names(authors.table))
authors.counts<-sapply(1:length(authors.match), function(x) ifelse(is.na(authors.match[x]),0,authors.table[authors.match[x]]))

# Chi-square test for goodness of fit to Poisson
authors.fit<-goodfit(authors.degree,type="poisson", par=list(lambda=mean(authors.degree)))

# Create df
authors.dd<-data.frame(list(val=0:max(authors.degree), counts=authors.counts, dd=authors.dd, fitted=authors.fit$fitted))

authors.ddplot<-ggplot(authors.dd, aes(x=val))+geom_line(aes(y=dd, colour="Degree Distribution"))+
    geom_point(aes(y=dd,colour="Degree Distribution"))+xlab("Degree")+ylab("Frequency")+theme_bw()+
    scale_colour_manual(values=c("Degree Distribution"="black"),name="Author Network")+
    opts(legend.position=c(.85,.85),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(plot=authors.ddplot,filename="../images/authors_ddplot.pdf",height=5,width=8)

authors.fitplot<-ggplot(authors.dd, aes(x=val))+geom_rect(aes(xmin=val-0.5,xmax=val+0.5,ymin=fitted-counts,ymax=fitted,fill="lightgrey",colour="Observed"))+
    geom_line(aes(y=fitted, colour="Fitted Poisson"))+geom_point(aes(y=fitted, colour="Fitted Poisson"))+
    geom_hline(yintercept=0,linetype=2,alpha=0.5)+theme_bw()+scale_colour_manual(values=c("Fitted Poisson"="black","Observed"="darkgrey"),name="Author Network")+
    scale_fill_manual(values=c("lightgrey"="lightgrey"),legend=FALSE)+xlab("Degree")+ylab("")+scale_alpha(legend=FALSE)+
    opts(legend.position=c(.85,.85),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(plot=authors.fitplot,filename="../images/authors_fitplot.pdf",height=5,width=8)

# Article network degree distribution
articles.degree<-degree(articles)
articles.dd<-degree.distribution(articles)
articles.table<-table(articles.degree)
articles.match<-match(0:max(articles.degree),row.names(articles.table))
articles.counts<-sapply(1:length(articles.match), function(x) ifelse(is.na(articles.match[x]),0,articles.table[articles.match[x]]))

# Chi-square test for goodness of fit to Poisson
articles.fit<-goodfit(articles.degree,type="poisson", par=list(lambda=mean(articles.degree)))

# Create data frame and plot
articles.dd<-data.frame(list(val=0:max(articles.degree), counts=articles.counts, dd=articles.dd, fitted=articles.fit$fitted))

articles.ddplot<-ggplot(articles.dd, aes(x=val))+geom_line(aes(y=dd, colour="Degree Distribution"))+
    geom_point(aes(y=dd,colour="Degree Distribution"))+xlab("Degree")+ylab("Frequency")+theme_bw()+
    scale_colour_manual(values=c("Degree Distribution"="black"),name="Article Network")+
    opts(legend.position=c(.85,.85),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(plot=articles.ddplot,filename="../images/articles_ddplot.pdf",height=5,width=8)

articles.fitplot<-ggplot(articles.dd, aes(x=val))+geom_rect(aes(xmin=val-0.5,xmax=val+0.5,ymin=fitted-counts,ymax=fitted,fill="lightgrey",colour="Observed"))+
    geom_line(aes(y=fitted, colour="Fitted Poisson",alpha=0.5))+geom_point(aes(y=fitted, colour="Fitted Poisson"))+
    geom_hline(yintercept=0,linetype=2)+theme_bw()+scale_colour_manual(values=c("Fitted Poisson"="black","Observed"="darkgrey"),name="Article Network")+
    scale_fill_manual(values=c("lightgrey"="lightgrey"),legend=FALSE)+xlab("Degree")+ylab("")+scale_alpha(legend=FALSE)+
    opts(legend.position=c(.85,.85),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(plot=articles.fitplot,filename="../images/articles_fitplot.pdf",height=5,width=8)

# Display results from Chi-squared test for 
author.fittest<-chisq.test(authors.dd$counts,authors.dd$fitted)
print(author.fittest)

article.fittest<-chisq.test(articles.dd$counts,articles.dd$fitted)
print(article.fittest)

# Both cases we cannot reject the null hypothesis that both distribution
# were generated with a Poisson, but the p-value provides evidence
# that the fit is at best weak.