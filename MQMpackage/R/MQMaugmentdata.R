#####################################################################
#
# MQMaugmentdata.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: MQMaugment
#
######################################################################

######################################################################
#
# scanMQM:
#
######################################################################
setwd("D:/")
library(qtl)
dyn.load("scanMQM.dll")
qtl <- c(3,15,3,7)							# QTL at chromosome 3
data(map10)									# Mouse genome
cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross

MQMaugment <- function(cross= NULL,maxaug=10000,maxiaug=1000,neglect=1000){
    library(qtl)
	if(is.null(cross)){
		print(paste("Error: No cross file. Please supply a valid cross object."))
		return 
	}
	if(class(cross)[1] == "f2"){
		#print(paste("INFO: Received an F2 cross."))
		n.ind <- nind(cross)
		n.chr <- nchr(cross)
		n.aug <- maxaug
		#print(paste("INFO: Number of individuals:",n.ind))
		#print(paste("INFO: Number of chr:",n.chr))
		geno <- NULL
		chr <- NULL
		dist <- NULL
		out.qtl <- NULL	
		for(i in 1:n.chr) {
			geno <- cbind(geno,cross$geno[[i]]$data)
			chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
			dist <- c(dist,cross$geno[[i]]$map)
		}
		pheno <- cross$pheno
		n.mark <- ncol(geno)
		#print(paste("INFO: Number of markers:",n.mark))
		for(i in 1:n.ind) {
			for(j in 1:n.mark) {
				if(is.na(geno[i,j])){
					geno[i,j] <- 9;
				}
			}
		}
		#check for missing phenotypes
		dropped <- NULL
		for(i in 1:dim(pheno)[1]) {
			if(is.na(pheno[i,1])){
			  print(paste("Dropped individual ",i," with missing genotype\n",sep=""))
			  dropped <- c(dropped,i) 
			  n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped,]
		}

		result <- .C("R_augdata",
                as.integer(geno),
				as.double(dist),
				as.double(pheno[,1]),
				augGeno=as.integer(rep(0,n.mark*maxaug)),
				augPheno=as.double(rep(0,maxaug)),
				as.integer(n.ind),
				as.integer(n.aug),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(maxaug),
				as.integer(maxiaug),
				as.integer(neglect)
			    )
		n.aug = result[[7]]
		#print(paste("INFO: Number of individuals:",n.aug))		
		markONchr <- 0
		markdone <- 0
		for(c in 1:n.chr){
			#print(paste("Cromosome",c,"\n",sep=""))
			matri <- NULL
			markONchr <- dim(cross$geno[[c]]$data)[2]
			#print(paste("# markers",markONchr,"\n",sep=""))
			for(j in markdone:(markdone+markONchr-1)){
			    #print(paste("Start",markdone,":End",(markdone+markONchr-1),"\n",sep=""))
				ind <- NULL
				pheno <- NULL
				for(i in 1:n.aug){
					ind = c(ind,result[[4]][i+(j*maxaug)])
					pheno <- rbind(pheno,result[[5]][i])
				}
				matri <- rbind(matri,ind)
			}
			matri <- t(matri)
			if(markdone==0){
			colnames(matri) <- colnames(geno)[markdone:(markdone+markONchr)]
			}else{
			#print(paste("Markdone",markdone,"End",(markdone+markONchr-1)))
			colnames(matri) <- colnames(geno)[(markdone+1):(markdone+markONchr)]			
			}
			cross$geno[[c]]$data <- matri
			colnames(pheno) = "phenotype"
			cross$pheno <- as.data.frame(pheno)

			markdone <- (markdone+markONchr)  
		}
		cross
	}else{
		print(paste("Error: Currently only f2 crosses can be analyzed by MQM."))
		return 
	}			
}

# end of MQMaugment.R

res <- MQMaugment(cross,maxaug=500,maxiaug=10,neglect=10)
