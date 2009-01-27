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

scanMQM <- function(cross){
    setwd("D:/")
    library(qtl)
	dyn.load("scanMQM.dll")
    cross <- sim.cross(map10,c(3,15,7,8),n=10)
    n.chr <- nchr(cross)
	geno <- NULL
	chr <- NULL
	for(i in 1:n.chr) {
      geno <- cbind(geno,cross$geno[[i]]$data)
	  chr <- c(chr,rep(i,dim(cross$geno[[1]]$data)[2]))
	}
	pheno <- cross$pheno
	n.ind <- nind(cross)
	n.mark <- ncol(geno)
	n.ind
	n.mark
	geno[1:5,1:5]
	result <- .C("R_scanMQM",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(1),    # 1 family
                as.integer(geno),
				as.integer(chr),
				as.double(pheno[,])
			    )
				
}

# end of scanMQM.R