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
dyn.load("scanMQM.dll")
cross <- read.cross("csv","","Test.csv")

scanMQM <- function(cross= NULL,cofactors = NULL,Phenot=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.run=0){
    library(qtl)
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2"){
		cat("INFO: Received an F2 cross.\n")
		n.ind <- nind(cross)
		n.chr <- nchr(cross)
		cat("INFO: Number of individuals: ",n.ind,"\n")
		cat("INFO: Number of chr: ",n.chr,"\n")
		geno <- NULL
		chr <- NULL
		dist <- NULL
		out.qtl <- NULL	
		for(i in 1:n.chr) {
			geno <- cbind(geno,cross$geno[[i]]$data)
			chr <- c(chr,rep(i,dim(cross$geno[[i]]$data)[2]))
			dist <- c(dist,cross$geno[[i]]$map)
		}
		if(alfa <=0 || alfa >= 1){
			stop("Error: Alfa must be between 0 and 1.")
		}
		if(n.run <0 || n.run >= 10000){
			stop("Error: # of runs should be positive and < 10000.")
		}
		pheno <- cross$pheno
		n.mark <- ncol(geno)
		cat("INFO: Number of markers:",n.mark,"\n")
		for(i in 1:n.ind) {
			for(j in 1:n.mark) {
				if(is.na(geno[i,j])){
					stop("ERROR: Missing genotype information, please estimate unknown data, before running scanMQM.\n")
					geno[i,j] <- 9
				}
			}
		}
		#check for missing phenotypes
		dropped <- NULL
		for(i in 1:dim(pheno)[1]) {
			if(is.na(pheno[i,1])){
			  cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			  dropped <- c(dropped,i) 
			  n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped,]
		}
		if(!is.null(cross$extra)){
			cat("INFO: previously augmented dataset.\n")			
			cat("INFO: Individuals before augmentation",cross$extra$Nind,".\n")
			extra1 <- cross$extra$Nind
			cat("INFO: Individuals after augmentation",cross$extra$augIND,".\n")
			extra2 <- cross$extra$augIND
		}else{
			extra1 <- n.ind
			extra2 <- 0:n.ind
		}
		#check if we have cofactors, so we can do backward elimination
		backward <- 0;
		if(is.null(cofactors)){
			cat("INFO: No cofactors, setting cofactors to 0\n")
			cofactors = rep(0,n.mark)
		}else{
			if(length(cofactors) != n.mark){
				cat("Error: # Cofactors != # Markers\n")		
			}else{
				print(paste("INFO:#",length(cofactors),"Cofactors received"),quote = FALSE)
				if(sum(cofactors) > 0){
					cat("INFO: Doing backward elimination of selected cofactors.\n")
					cat("INFO: Not doing permutation of data.\n")
					backward <- 1;
					n.run <- 0;
				}else{
					backward <- 0;
					stop("Error: Are u trying to give an empty cofactor list ???")
				}
			}
		}
		if(Phenot != 1){
			cat("INFO: Selected phenotype ",Phenot,".\n")
			cat("INFO: # of phenotypes in object ",nphe(cross),".\n")
			if(nphe(cross) < Phenot || Phenot < 1){
				stop("Error: No such phenotype")
			}			
		}
		qtlAchromo <- length(seq(step.min,step.max,step.size))
		cat("Number of locations per chromosome: ",qtlAchromo, "\n")
		result <- .C("R_scanMQM",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
                as.integer(geno),
				as.integer(chr),
				as.double(dist),
				as.double(pheno[,Phenot]),
				as.integer(cofactors),
				as.integer(backward),
				as.integer(REMLorML),
				as.double(alfa),
				as.integer(em.iter),
				as.double(windowsize),
				as.double(step.size),
				as.double(step.min),
				as.double(step.max),
				as.integer(n.run),
				as.integer(extra1),
				as.integer(extra2),
				QTL=as.double(rep(0,n.chr*qtlAchromo))
			    )
		# initialize output object
		qtl <- NULL
		names <- NULL
		for(i in 1:(n.chr*qtlAchromo)) {
			qtl <- rbind(qtl,c(ceiling(i/qtlAchromo),rep(seq(step.min,step.max,step.size),n.chr)[i],result$QTL[i]))
			names <- c(names,paste("C",ceiling(i/qtlAchromo),"L",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
		}
		rownames(qtl) <- names
		colnames(qtl) = c("chr","pos (Cm)","QTL")
		
		#So we can use carls plotting routines
		class(qtl) <- c(class(qtl),"scanone") 
	
		#return QTL
		qtl
	
	}else{
		stop("Error: Currently only f2 crosses can be analyzed by MQM.")
	}			
}

# end of scanMQM.R

res <- scanMQM(cross)
plot(res)
