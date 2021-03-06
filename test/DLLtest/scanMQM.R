#####################################################################
#
# scanMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: scanMQM
#
######################################################################

######################################################################
#
# scanMQM:
#
######################################################################
setwd("D:/")
library(qtl)
cross <- read.cross("csv","","Test.csv")

scanMQM <- function(cross= NULL,cofactors = NULL,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0){
    library(qtl)
	dyn.load("scanMQM.dll")
	if(is.null(cross)){
      print(paste("Error: No cross file. Please supply a valid cross object."))
	  return 
	}
if(class(cross)[1] == "f2"){
    print(paste("INFO: Received an F2 cross."))
	n.ind <- nind(cross)
	n.chr <- nchr(cross)
	print(paste("INFO: Number of individuals:",n.ind))
	print(paste("INFO: Number of chr:",n.chr))
	geno <- NULL
	chr <- NULL
	dist <- NULL
	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	  chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
	  dist <- c(dist,cross$geno[[i]]$map)
	}
	if(alfa <=0 || alfa >= 1){
      print(paste("Error: Alfa must be between 0 and 1."))
	  return	
	}
	pheno <- cross$pheno
	n.mark <- ncol(geno)
	print(paste("INFO: Number of markers:",n.mark))
	for(i in 1:n.ind) {
	  for(j in 1:n.mark) {
	    if(is.na(geno[i,j])){
	      geno[i,j] <- 9;
	    }
	  }
	}
	backward <- 0;
	if(is.null(cofactors)){
	 print(paste("INFO: No cofactors, setting cofactors to 0"))
	 cofactors = rep(0,n.mark)
	}else{
	  if(length(cofactors) != n.mark){
	    print("Error: # Cofactors != # Markers")		
	  }else{
	    print(paste("INFO:#",length(cofactors),"Cofactors received"))
		if(sum(cofactors) > 0){
		  print(paste("INFO: Doing backward elimination of selected cofactors."))
		  backward <- 1;
		}else{
		  backward <- 0;
		  print(paste("Error: Are u trying to give an empty cofactor list ???"))
		}
	  }
	}

	result <- .C("R_scanMQM",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(1),    # 1 family
                as.integer(geno),
				as.integer(chr),
				as.double(dist),
				as.double(pheno[,1]),
				as.integer(cofactors),
				as.integer(backward),
				as.integer(REMLorML),
				as.double(alfa),
				as.integer(em.iter),
				as.double(windowsize),
				as.double(step.size),
				as.double(step.min),
				as.double(step.max)
			    )
}else{
    print(paste("Error: Currently only f2 crosses can be analyzed by MQM."))
	return 
}
				
}

# end of scanMQM.R

scanMQM(cross)
