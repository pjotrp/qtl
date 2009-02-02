#####################################################################
#
# makeCofList.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: makeCofList
#
######################################################################

######################################################################
#
# makeCofList: Prepares a cofactor list to use with scanMQM
#
######################################################################

setwd("D:/")
library(qtl)
dyn.load("scanMQM.dll")
cross <- read.cross("csv","","Test.csv")

makeCofList <- function(cross= NULL,cofactors = NULL,sexfactors=NULL){
	if(is.null(cross)){
      print(paste("Error: No cross file. Please supply a valid cross object."))
	  return 
	}
	if(is.null(cofactors)){
      print(paste("Error: Cofactors to set. Please supply a list of markers to serve as cofactors."))
	  return 
	}
		
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL

	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	
	n.mark <- ncol(geno)

	if(max(cofactors) > n.mark){
      print(paste("Error: Trying to set a non-existent marker as a cofactor."))
	  return 	  
	}

	if(max(sexfactors) > n.mark){
      print(paste("Error: Trying to set a non-existent marker as a sexfactor."))
	  return 	  
	}

    cofactorlist <- rep(0,n.mark)
	for(i in 1:length(cofactors)) {
	  cofactorlist[cofactors[i]]=1
	}
	if(!is.null(sexfactors)){
      for(i in 1:length(sexfactors)) {
	    cofactorlist[sexfactors[i]]=2
	  }
	}
    cofactorlist
}

a <- makeCofList(cross,c(10,20,30,40,50,60,70,80),c(186,187))
b <- makeCofList(cross,c(10,20,30,40,50,60,70,80,300),c(186,187))