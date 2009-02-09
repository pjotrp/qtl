#####################################################################
#
# MQMCofactors.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: MQMCofactors
#
######################################################################

######################################################################
#
# MQMCofactors: Prepares a cofactor list to use with scanMQM
#
######################################################################

setwd("D:/")
library(qtl)
dyn.load("scanMQM.dll")
cross <- read.cross("csv","","Test.csv")

MQMCofactors <- function(cross= NULL,cofactors = NULL,sexfactors=NULL){
	if(is.null(cross)){
      stop("Error: No cross file. Please supply a valid cross object.")
	  return 
	}
	if(is.null(cofactors)){
      stop("Error: Cofactors to set. Please supply a list of markers to serve as cofactors.")
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
      stop("Error: Trying to set a non-existent marker as a cofactor.")
	  return 	  
	}

	if(min(cofactors) <= 0){
      stop("Error: Trying to set a non-existent marker as a cofactor.")
	  return 	  
	}
	
	if(!is.null(sexfactors)){
	if(max(sexfactors) > n.mark){
      stop("Error: Trying to set a non-existent marker as a sexfactor.")
	  return 	  
	}
	
	if(min(sexfactors) <= 0){
      stop("Error: Trying to set a non-existent marker as a sexfactor.")
	  return 	  
	}
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

MQMCofactorsEach <- function(cross = NULL,each = 3){
	if(is.null(cross)){
      stop("Error: No cross file. Please supply a valid cross object.")
	  return 
	}


	if(each < 2){
      stop("Error: Can't set cofactors that often.")
	  return 
	}

	
	n.chr <- nchr(cross)
	geno <- NULL
	cofactorlist <- NULL

	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	}
	
	n.mark <- ncol(geno)

	if(each > n.mark){
      stop("Error: Not enough markers to place cofactors at.")
	  return 
	}	
	
    cofactorlist <- rep(0,n.mark)
	for(i in 1:n.mark) {
		if(i%%as.integer(each)==0){
			cofactorlist[i] = 1
		}
	}
    cofactorlist
}

a <- MQMCofactors(cross,c(10,20,30,40,50,60,70,80),c(186,187))