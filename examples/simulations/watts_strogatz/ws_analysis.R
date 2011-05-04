# # File-Name:       ws_analysis.R           
# # Date:            2011-04-11                                
# # Author:          Drew Conway
# # Email:           drew.conway@nyu.edu                                      
# # Purpose:         Statistically explore the results of the Watts-Strogatz GMM models
# # Data Used:       Output from watts_strogatz.py
# # Packages Used:   igraph, ggplot2
# # Output File:    
# # Data Output:     
# # Machine:         Drew Conway's MacBook Pro
# 
# # Copyright (c) 2011, under the Simplified BSD License.  
# # For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
# # All rights reserved.                                                         
# 
library(ggplot2)
library(igraph)

# Base WS model data
ws.df<-read.csv("ws_df.csv", stringsAsFactors=FALSE)

# Normalizing constants
L0<-50.450450450450454
C0<-0.6666666666666636
probs.breaks<-c(.0001, .001, .01, .1, 1)

# Take means of C(p) and L(p) and normalize
ws.norm<-ddply(ws.df, .(p), summarise, C.norm=mean(cc)/C0, L.norm=mean(avg.pl)/L0)
    
# Simulated GMM WS model data
sims.df<-read.csv("sims_df.csv", stringsAsFactors=FALSE)

# Take means of C(p) and L(p) and normalize
sims.norm<-ddply(sims.df, .(p), summarise, C.norm=mean(cc)/C0, L.norm=mean(avg.pl)/L0)

# Create plot to compare the two mocels using clustering coeffecient, C(p), and 
# average shortest path length, L(p)
ws.comp<-ggplot(ws.norm, aes(x=p))+geom_point(aes(y=C.norm, color="C(p)/C(0)"))+
    geom_point(aes(y=L.norm, color="L(p)/L(0)"))+
    geom_line(aes(y=C.norm, color="C(p)/C(0)", linetype="WS"))+
    geom_line(aes(y=L.norm, color="L(p)/L(0)",linetype="WS"))+
    geom_point(data=sims.norm, aes(y=C.norm, color="C(p)/C(0)"))+
    geom_point(data=sims.norm, aes(y=L.norm, color="L(p)/L(0)"))+
    geom_line(data=sims.norm, aes(y=C.norm, color="C(p)/C(0)", linetype="GMM"))+
    geom_line(data=sims.norm, aes(y=L.norm, color="L(p)/L(0)",linetype="GMM"))+
    scale_color_manual(values=c("C(p)/C(0)"="darkred","L(p)/L(0)"="darkblue"), name="Statistic")+
    scale_linetype_manual(values=c("WS"=1,"GMM"=2), name="Model type")+
    scale_x_log10(breaks=probs.breaks, labels=probs.breaks)+
    scale_y_continuous(breaks=seq(0,1,.1))+
    opts(title="Comparison Between Watts-Strogatz Random Graph\nModel and GMM Equivalent [n=100, k=3]")+
    xlab("p")+ylab("")+
    theme_bw()
    
# Save output
ggsave(plot=ws.comp, filename="ws_comp.pdf", height=7, width=10)