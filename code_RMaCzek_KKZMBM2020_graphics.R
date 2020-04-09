## Copyright (C) 2020 
## Author Krzysztof Bartoszek
##
##  This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program. If not, see <https://www.gnu.org/licenses/>.

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## script generating the graphics for the manuscript
## the axes labels in the manuscript are manually edited (rewritten in larger fonts) for visibility


load("results_RMaCzek_KKZMBM2020_JCze1909setup.RData")
## compare with Czekanowski (1909)'s original distances
mDistJCze<-as.matrix(res_all[[1]]$mDistJCzeTraitsSkulls_distf_direct)
skulls_distances_2<-skulls_distances ## Czekanowski (1909)'s distance matrix is part of RMaCzek
skulls_distances_2["Galey Hill","Neandertal"]<-10.540 ## here there was a typo in Czekanowski (1909)'s data
rownames(skulls_distances_2)[3]<-"Krapina C" ## here there was a typo in Czekanowski (1909)'s data
dsdDJ<-(abs(skulls_distances_2-mDistJCze))/skulls_distances_2 ## the relative difference between the two distance matrices
diag(dsdDJ)<-0 ## diagonal is 0
png("plot_RMaCzek_KKZMBM2020_htmpdsdDJ.png");heatmap.2(dsdDJ,Rowv=NA,Colv=NA,scale="none", dendrogram="none",density.info = "none",trace = "none",lhei = c(12,2),sepwidth =c(0.001,0.0001),symkey=FALSE,key.title=NA,key.xlab=NA,lmat=rbind( c(1,3),c(2,4)),lwid=c(7,2),col=terrain.colors);dev.off()
## plot Czekanowski's diagram under the best found ordering for the replication of Czekanowski (1909)'s analysis
png("plot_RMaCzek_KKZMBM2020_JCze1909ordered.png");plot(res_all[[1]]$res_distf_from_direct,plot_title="");dev.off() 


## plot some of Czekanowski's diagrams for the other setups
load("results_RMaCzek_KKZMBM2020_OtherSetups.RData")
## draw example figures for illustrating the results
v_indextodraw<-c(13,21,25,26,29,46,66,82,92) ## example ones to plot, corresponding to specific placements of the Nowsiolka skull, bar nr. 46 which is interesting with respect to the Brux skull

for (i in v_indextodraw){
    png(paste0("plot_RMaCzek_KKZMBM2020_figsetup_",i,".png"))
    plot(res_all[[i]]$res_distf_from_direct,plot_title=paste0("setup ",i),cex.main=2)
    dev.off()
}

