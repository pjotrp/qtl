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

scanMQM <- function(cross= NULL,cofactors = NULL){
    library(qtl)
	dyn.load("scanMQM.dll")
	if(is.null(cross)){
      print(paste("No cross, generating a random cross from map10 dataset (100 individuals)"))
	  cross <- sim.cross(map10,n=100)
	}
	n.ind <- nind(cross)
	n.chr <- nchr(cross)
	print(paste("Number of individuals:",n.ind))
	print(paste("Number of chr:",n.chr))
	geno <- NULL
	chr <- NULL
	dist <- NULL
	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	  chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
	  dist <- c(dist,cross$geno[[i]]$map)
	}
	
	pheno <- cross$pheno
	n.mark <- ncol(geno)
	print(paste("Number of markers:",n.mark))
	backward <- 0;
		
	
	if(is.null(cofactors)){
	 print(paste("No cofactors, setting cofactors to 0"))
	 cofactors = rep(0,n.mark)
	}else{
	  if(length(cofactors) != n.mark){
	    print("Error: # Cofactors != # Markers")		
	  }else{
	    print(paste("Cofactors received"))
		if(sum(cofactors) > 0){
		  backward <- 1;
		}else{
		  backward <- 0;
		  print(paste("Sneaky, are u trying to give an empty cofactor list ???"))
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
				as.integer(backward)
			    )
				
}

# end of scanMQM.R