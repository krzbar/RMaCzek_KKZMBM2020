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

## References
## [1] J. Czekanowski (1909). Zur Differentialdiagnose der Neandertalgruppe, Korespondentblatt der
## Deutschen Gesellschaft fur Anthropologie, Ethnologie und Urgeschichte, XL(6/7), 44-47.
## [2] K. Stolyhwa (1908). Czaszka z Nowosiolki jako dowod istnienia w okresie historycznym kształtow
## pokrewnych z Homo primigenius, Rozprawy Wydzialu matematyczno-przyrodnicznego
## Akademii Umiejetnosci, XLVIII(B), 1-27 (The skull from Nowosiolka  as proof of existence during the era of 
## history shapes common with Homo primigenius, in Polish).

library(RMaCzek) ## requires RMaCzek >= 1.3.3 to run
library(vegan)
library(gplots)

source("code_RMaCzek_KKZMBM2020_distancefunctions.R") ## file with functions calculating the distances

b_use_random_seed_from_manuscript<-TRUE ## set to FALSE if an independent run is wanted and not a replication of the results in the manuscript

rexp(1) ## initialize the a random number generator
for (i in 1:2){
    if (i==1){
	## Settings reanalyzing Czekanowski (1909)'s paper
	runnum<-"JCze1909setup_newrun"
	v_b_use_JCzeTrait<-c(TRUE) ## use only the 27 traits used by Czekanowski (1909)
	v_b_radians<-c(FALSE) ## degrees will not be changed into radians
	v_b_remove_ratio_variables<- c(0) ## do not remove index variables
	v_num_similar<-c(3) ## setting ignored
	v_b_use_mean_Sol<-c(TRUE)  ## setting ignored
	v_c_dist_func<-c("DD") ## use Czekanowski (1909)'s average distance function
	v_NA_frac_obs<-c(NA) ## setting ignored
	v_NA_frac_var<-c(NA) ## setting ignored
	b_do_scale<-FALSE ## do not normalize, mean center, divide by standard deviation, the data, as not done by Czekanowski (1909)
	if(b_use_random_seed_from_manuscript){## load the random seed to replicate the result
	    load("randomseed_RMaCzek_KKZMBM2020_JCze1909setup.RData") 
	}
    }

    if (i==2){
	## All settings to try out different possibilities
	runnum<-"OtherSetups_newrun"
	v_b_use_JCzeTrait<-c(TRUE,FALSE) ## use only the 27 traits (TRUE) used by Czekanowski (1909) or all (FALSE)
	v_b_radians<-c(TRUE,FALSE) ## degrees changed into radians (TRUE) or not (FALSE)
	v_b_remove_ratio_variables<- c(0,1,2) ##0: all data as is, 1: remove variables that are a ration of two other, 2: keep the ratio variables but remove its two components
	v_num_similar<-c(3) ## for Stolyhwa distance: should a given variable be classfied to which of two representative observations it is closer to (2 levels) or also in-between (3 levels)
	v_b_use_mean_Sol<-c(TRUE,FALSE) ## for Stolyhwa distance: should the representative observation be an average of the observations in the class
	v_c_dist_func<-c("Stolyhwa","DD", "L2") ## distance function. Stolyhwa: counts the number of traits that differntiate between to observations when comparing to representative observations from each class, divided by number of variables; DD: average difference Czekanowski (1909); L2: squared Euclidean distance
	v_NA_frac_obs<-c(NA,0.5) ##all observations taken, or those with missingness on more than 50% variables removed
	v_NA_frac_var<-c(NA) ## all variables considered, here one can p
	b_do_scale<-TRUE ## normalize, mean center, divide by standard deviation, the data
	if(b_use_random_seed_from_manuscript){## load the random seed to replicate the results
	    load("randomseed_RMaCzek_KKZMBM2020_OtherSetups.RData") ## load the random seed to replicate the results
	}
    }

    if(b_use_random_seed_from_manuscript){##setup the random number generator if replicating results
	RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
	RNGversion(RNG_version)     
    }

    ## all possible setups
    df_all_setups<-expand.grid(b_use_JCzeTrait=v_b_use_JCzeTrait,b_radians=v_b_radians,b_remove_ratio_variables=v_b_remove_ratio_variables,num_similar=v_num_similar,b_use_mean_Sol=v_b_use_mean_Sol,c_dist_func=v_c_dist_func,NA_frac_obs=v_NA_frac_obs,NA_frac_var=v_NA_frac_var)

    ## remove setups that are the same, i.e. some control variables are ignored by the L2 and DD distance functions
    v_torem<-intersect(which(df_all_setups[,"c_dist_func"]=="DD"),which(df_all_setups[,"num_similar"]==2))
    v_torem<-c(v_torem,intersect(which(df_all_setups[,"c_dist_func"]=="DD"),which(df_all_setups[,"b_use_mean_Sol"]==FALSE)))
    v_torem<-c(v_torem,intersect(which(df_all_setups[,"c_dist_func"]=="L2"),which(df_all_setups[,"num_similar"]==2)))
    v_torem<-c(v_torem,intersect(which(df_all_setups[,"c_dist_func"]=="L2"),which(df_all_setups[,"b_use_mean_Sol"]==FALSE)))
    if (length(v_torem)>0){df_all_setups<-df_all_setups[-v_torem,,drop=FALSE]}

    c_rep_HS<-"Brux" ## H. sapiens skull, has all 27, used by Czekanowski (1909) variables observed, taken as representative skull for Stolyhwa disance
    c_rep_Nead<-"Neandertal" ## H. neandertalsis skull, has all 27, used by Czekanowski (1909) variables observed, taken as representative skull for Stolyhwa disance
    v_HS<-NA

    ## Nowosiolka is not used as it is meant to be the focal species that we want to place
    v_HS<-c("Kannstatt", "Galey.Hill", "Brunn", "Brux", "Egisheim") ## H. sapiens skulls for average skull 
    v_Nead<-c("Spy.I", "Spy.II", "Krapina.C", "Krapina.D", "Neandertal", "Gibraltar", "Pithecanthropus") ## H. neandertalsis skulls for average skull 

    ## Read in the data presented by Stolyhwa (1908)
    dfHumanoidData<-read.csv("data_RMaCzek_KKZMBM2020_Stolyhwa1908_craniometricdata.csv",header=TRUE,sep=";")


    ## Create data frame for further analyses
    v_skullnames<-colnames(dfHumanoidData)[4:33]
    v_traitnames<-c()
    v_JCzIndex<-c()
    v_kept_row_ids<-c()
    dfavgHumanoidData<-matrix(NA,ncol=length(v_skullnames),nrow=1)
    j<-1
    ## Stolyhwa (1908) reports some variables as single number but others as a range (min,max)
    ## here we identify which ones are reported as a single number (copied directly)
    ## and which as a range, then the average (min+max)/2 is taken
    for (i in 1:(nrow(dfHumanoidData))){
	if ((i<nrow(dfHumanoidData))&&(dfHumanoidData[i,1]=="trait_value")&&(dfHumanoidData[i+1,1]=="trait_min")){
	    v_vals<-c(as.matrix(dfHumanoidData[i,v_skullnames]))
	    v_isNA<-which(is.na(v_vals))
	    v_vals[v_isNA]<-((c(as.matrix(dfHumanoidData[i+1,v_skullnames]))+c(as.matrix(dfHumanoidData[i+2,v_skullnames])))/2)[v_isNA]
	    dfavgHumanoidData<-rbind(dfavgHumanoidData,rep(NA,length(v_skullnames)))
	    dfavgHumanoidData[j,]<-v_vals
	    v_traitnames<-c(v_traitnames,levels(dfHumanoidData[i,"trait_name_pl"])[dfHumanoidData[i,"trait_name_pl"]])
	    v_JCzIndex<-c(v_JCzIndex,dfHumanoidData[i,"Czekanowskis_trait_index_TabelleI"])
	    j<-j+1
	    v_kept_row_ids<-c(v_kept_row_ids,i)	
	}else{
	    if (((dfHumanoidData[i,]=="trait_value")&&(dfHumanoidData[i+1,1]=="trait_value"))){
		dfavgHumanoidData<-rbind(dfavgHumanoidData,rep(NA,length(v_skullnames)))
		dfavgHumanoidData[j,]<-c(as.matrix(dfHumanoidData[i,v_skullnames]))
		v_traitnames<-c(v_traitnames,levels(dfHumanoidData[i,"trait_name_pl"])[dfHumanoidData[i,"trait_name_pl"]])
		v_JCzIndex<-c(v_JCzIndex,dfHumanoidData[i,"Czekanowskis_trait_index_TabelleI"])
		j<-j+1
		v_kept_row_ids<-c(v_kept_row_ids,i)	
	    }
	}
    }

    m_trait_names<-cbind(names_index=1:length(v_kept_row_ids),dfHumanoidData[v_kept_row_ids,c("units","Czekanowskis_trait_index_TabelleI","components_index_variable","dependent_index_variable"),drop=FALSE],b_deg=FALSE)
    m_trait_names[which(m_trait_names[,"units"]=="degree"),"b_deg"]<-TRUE

    dfavgHumanoidData<-dfavgHumanoidData[-nrow(dfavgHumanoidData),,drop=FALSE] ## remove last row as a new one is always added to the end in the previous for loop
    colnames(dfavgHumanoidData)<-v_skullnames
    rownames(dfavgHumanoidData)<-v_traitnames
    dfavgHumanoidData<-cbind(dfavgHumanoidData,JCze1909index=v_JCzIndex) ## vector telling us if the variables was used by Czekanowski (1909)

    res_all<-apply(df_all_setups,1,function(v_setup,dfavgHumanoidData,m_trait_names,b_do_scale){
    ## run the analysis for each setup
	b_use_JCzeTrait<-v_setup["b_use_JCzeTrait"]
	b_radians<-v_setup["b_radians"]
        b_remove_ratio_variables<-v_setup["b_remove_ratio_variables"]
	num_similar<-v_setup["num_similar"]
	b_use_mean_Sol<-v_setup["b_use_mean_Sol"]
	c_dist_func<-v_setup["c_dist_func"]
	NA_frac_var<-as.numeric(v_setup["NA_frac_var"])
	NA_frac_obs<-as.numeric(v_setup["NA_frac_obs"])

	if (b_remove_ratio_variables==1){
	## we do not want Egisheim as it for some has only ration measrurements, missing the raw ones
	    v_HS<-c("Kannstatt", "Galey.Hill", "Brunn", "Brux")
	}

	if (b_radians){
    	    dfavgHumanoidData[which(m_trait_names[,"b_deg"]),]<-2*pi*dfavgHumanoidData[which(m_trait_names[,"b_deg"]),]/360
	}
	if (b_remove_ratio_variables==1){
	    dfavgHumanoidData<-dfavgHumanoidData[-which(m_trait_names[,"dependent_index_variable"]),]
	}
        if (b_remove_ratio_variables==2){
		dfavgHumanoidData<-dfavgHumanoidData[-which(m_trait_names[,"components_index_variable"]),]
	}

	dfavgHumanoidDataJCzeTraits<-dfavgHumanoidData
	if (b_use_JCzeTrait){
	    dfavgHumanoidDataJCzeTraits<-dfavgHumanoidData[which(!is.na(dfavgHumanoidData[,"JCze1909index"])),]
	}

	if (b_remove_ratio_variables==1){
	    ## we do not want Egisheim as it for some has only ration measrurements, missing the raw ones
	    v_JCzeSkulls<-c("Spy.I", "Spy.II", "Krapina.C", "Krapina.D", "Neandertal", "Gibraltar", "Pithecanthropus", "Kannstatt", "Galey.Hill", "Brunn", "Brux", "Nowosiolka")
	}else{
	    v_JCzeSkulls<-c("Spy.I", "Spy.II", "Krapina.C", "Krapina.D", "Neandertal", "Gibraltar", "Pithecanthropus", "Kannstatt", "Galey.Hill", "Brunn", "Brux", "Egisheim", "Nowosiolka")
	}
    
	## Keep only the skulls used by Czekanowski (1909), the others presented by Stolyhwa (1908) have nearly everything missing
	dfavgHumanoidDataJCzeTraitsSkulls<-dfavgHumanoidDataJCzeTraits[,v_JCzeSkulls]

        mDistJCzeTraitsSkulls_distf<-NA
	if (c_dist_func=="Stolyhwa"){## similarity by number of variables which are same or different	
	    f_dist_fun<-f_Sol1908dist
	    ## with this distance we cannot center and scale the data as the similarity to the representative observation could be lost
	    b_do_scale<-FALSE

	    if (b_use_mean_Sol){## average of observations used
		v1<-apply(dfavgHumanoidDataJCzeTraitsSkulls[,v_HS],1,mean,na.rm=TRUE)
		v2<-apply(dfavgHumanoidDataJCzeTraitsSkulls[,v_Nead],1,mean,na.rm=TRUE)
	    }else{## representative observation used
		v1<-dfavgHumanoidDataJCzeTraitsSkulls[,c_rep_HS]
		v2<-dfavgHumanoidDataJCzeTraitsSkulls[,c_rep_Nead]
	    }
	    ## list of arguments for the distance function to be passed to RMaCzek::czek_matrix()
	    dist_args<-list(dist_fun=f_Sol1908dist,NAfrac_obs=NA_frac_obs,NAfrac_var=NA_frac_var,b_do_scale=FALSE,v1=v1,v2=v2,num_similar=num_similar)	
	}

	if (c_dist_func=="DD"){## replication attempt at Czekanowski (1909)'s analysis
	    f_dist_fun<-f_JCze1909dist
	    v1<-NULL;v2<-NULL;num_similar<-NULL
	    dist_args<-list(dist_fun=f_JCze1909dist,NAfrac_obs=NA_frac_obs,NAfrac_var=NA_frac_var,b_do_scale=FALSE)
	}

	if (c_dist_func=="L2"){## use squared Euclidean distance
	    f_dist_fun<-f_L2
	    v1<-NULL;v2<-NULL;num_similar<-NULL
	    dist_args<-list(dist_fun=f_L2,NAfrac_obs=NA_frac_obs,NAfrac_var=NA_frac_var,b_do_scale=FALSE)
	}

	mDistJCzeTraitsSkulls_distf_direct<-NA
	mDistJCzeTraitsSkulls_distf_czek_matrix<-NA
	res_distf_from_direct<-NA
	res_distf_from_czek_matrix<-NA
	res_distf_from_data<-NA

	tryCatch({
	## Examples of various ways of obtaining the distances and diagrams
	    mDistJCzeTraitsSkulls_distf_direct<-f_calcdist(dfavgHumanoidDataJCzeTraitsSkulls,dist_fun=f_dist_fun,NAfrac_obs=NA_frac_obs,NAfrac_var=NA_frac_var,b_as_dist=TRUE,b_do_scale=b_do_scale,v1=v1,v2=v2,num_similar=num_similar)
	    mDistJCzeTraitsSkulls_distf_czek_matrix<-RMaCzek::czek_matrix(dfavgHumanoidDataJCzeTraitsSkulls,order=NA,distfun=f_calcdist,scale_data=b_do_scale,as_dist=TRUE,dist_args=dist_args)
	    res_distf_from_direct<-RMaCzek::czek_matrix(as.dist(mDistJCzeTraitsSkulls_distf_direct),order="OLO",scale_data=b_do_scale)
	    res_distf_from_czek_matrix<-RMaCzek::czek_matrix(as.dist(mDistJCzeTraitsSkulls_distf_czek_matrix),order="OLO",scale_data=b_do_scale)
	## below is the recommended way, where RMaCzek::czek_matrix does everything by itself
	    res_distf_from_data<-RMaCzek::czek_matrix(dfavgHumanoidDataJCzeTraitsSkulls,order="OLO",distfun=f_calcdist,scale_data=b_do_scale,dist_args=dist_args)
	},error=function(e){})
	res<-list(setup=v_setup,dfdata=dfavgHumanoidDataJCzeTraitsSkulls,mDistJCzeTraitsSkulls_distf_direct=mDistJCzeTraitsSkulls_distf_direct,mDistJCzeTraitsSkulls_distf_czek_matrix=mDistJCzeTraitsSkulls_distf_czek_matrix,res_distf_from_direct=res_distf_from_direct,res_distf_from_czek_matrix=res_distf_from_czek_matrix,res_distf_from_data=res_distf_from_data)

	res
    },dfavgHumanoidData=dfavgHumanoidData,m_trait_names=m_trait_names,b_do_scale=b_do_scale)

    save(res_all,file=paste0("results_RMaCzek_KKZMBM2020_",runnum,".RData"))

    ## Print to a text file the found ordering
    sink(paste0("results_RMaCzek_KKZMBM2020_",runnum,".txt"))
    for (i in 1:length(res_all)){
	cat(i);cat("\n")
	print(res_all[[i]]$setup)
	print("Data")
	print(res_all[[i]]$res_distf_from_data)
	print("=======================================")
    }
    sink()
}

## create the graphics for the manuscript (or equivalent ones from a new run)
source("code_RMaCzek_KKZMBM2020_graphics.R")
