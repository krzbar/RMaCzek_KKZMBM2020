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

## distance functions used for the skulls' analyses

f_JCze1909dist<-function(x,y,...){
## Czekanowski's average difference
    mean(abs(x-y),na.rm=TRUE)
}

f_L2<-function(x,y,...){
## squared Euclidean distance
    sum((x-y)^2,na.rm=TRUE)
}

f_Sol1908dist<-function(x,y,v1,v2,num_similar=2,...){
## distance function counts how many variables make a given observation closer to one representative observation (v1)
## how many to the other (v2) and (if num_similar==3) how many variables are in-between
## then the distance between the two observations is the mean of difference between the similarity counts of the variables

    v_use_names<-names(x)
    v_use_names<-intersect(v_use_names,names(y))
    v_use_names<-intersect(v_use_names,names(v1))
    v_use_names<-intersect(v_use_names,names(v2))
    v1<-v1[v_use_names];v2<-v2[v_use_names];x<-x[v_use_names];y<-y[v_use_names]
    
    x1<-abs(x-v1)
    x2<-abs(x-v2)
    y1<-abs(y-v1)
    y2<-abs(y-v2)
    if (num_similar==2){
	xd<-as.numeric(x1<x2)
	yd<-as.numeric(y1<y2)
    }
    if (num_similar==3){
	v_x<-x1-x2
	v_y<-y1-y2
	qx<-quantile(v_x,probs=c(1/3,2/3),na.rm=TRUE)
	qy<-quantile(v_x,probs=c(1/3,2/3),na.rm=TRUE)
	xd<-rep(NA,length(v_x))
	yd<-rep(NA,length(v_y))
	xd[which(v_x<=qx[1])]<-0
	xd[intersect(which(v_x>qx[1]),which(v_x<=qx[2]))]<-1
	xd[which(v_x>qx[2])]<-2

	yd[which(v_y<=qy[1])]<-0
	yd[intersect(which(v_y>qy[1]),which(v_y<=qy[2]))]<-1
	yd[which(v_y>qy[2])]<-2
	
    }

    res<-mean(abs(xd-yd),na.rm=TRUE)
    res
}


f_calcdist<-function(x,NAfrac_obs=NA,NAfrac_var=NA,dist_fun=f_JCze1909dist,b_as_dist=TRUE,b_do_scale=TRUE,...){
## a general function that will return the distance matrix under the user requested distance function
    if (b_do_scale){x<-scale(x)}
    if (!is.na(NAfrac_obs)){
	vremtr<-apply(x,1,function(x,NAfrac_obs){res<-FALSE;if(length(which(is.na(x)))>=(NAfrac_obs*length(x))){res<-TRUE};res},NAfrac_obs=NAfrac_obs)
	names(vremtr)<-NULL
	vremtr<-which(vremtr)
	if (length(vremtr)>0){x<-x[-vremtr,,drop=FALSE]}
    }
    if (!is.na(NAfrac_var)){
	vremtr<-apply(x,2,function(x,NAfrac_var){res<-FALSE;if(length(which(is.na(x)))>=(NAfrac_var*length(x))){res<-TRUE};res},NAfrac_var=NAfrac_var)
	names(vremtr)<-NULL
	vremtr<-which(vremtr)
	if (length(vremtr)>0){x<-x[,-vremtr,drop=FALSE]}
    }
    
    mDist<-matrix(NA,nrow=ncol(x),ncol=ncol(x))
    rownames(mDist)<-colnames(x)
    colnames(mDist)<-colnames(x)
    if (nrow(x)>0){
	for (i in 1:ncol(x)){
	    for (j in 1:ncol(x)){
		mDist[i,j]<-dist_fun(x[,i],x[,j],...)
	    }
	}
    }
    if (b_as_dist){mDist<-as.dist(mDist)}
    mDist
}

