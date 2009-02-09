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
data(listeria)

MQMaugment <- function(cross= NULL,Phenot=1,maxaug=1000,maxiaug=10,neglect=10){
    library(qtl)
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.")
		return 
	}
	if(class(cross)[1] == "f2"){
		
		n.ind <- nind(cross)
		n.chr <- nchr(cross)
		n.aug <- maxaug
		cat("INFO: Received an F2 cross.\n")
		cat("INFO: Number of individuals:",n.ind,"\n")
		cat("INFO: Number of chr:",n.chr,"\n")
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
		cat("INFO: Number of markers:",n.mark)
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
				cat("Dropped individual ",i ," with missing phenotype.\n")
				dropped <- c(dropped,i) 
				n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped,]
		}
		if(Phenot != 1){
			cat("INFO: Selected phenotype ",Phenot,".\n")
			cat("INFO: # of phenotypes in object ",nphe(cross),".\n")
			if(nphe(cross) < Phenot || Phenot < 1){
				stop("Error: No such phenotype")
			}			
		}
		result <- .C("R_augdata",
                as.integer(geno),
				as.double(dist),
				as.double(pheno[,Phenot]),
				augGeno=as.integer(rep(0,n.mark*maxaug)),
				augPheno=as.double(rep(0,maxaug)),
				augIND=as.integer(rep(0,maxiaug*n.ind)),
				as.integer(n.ind),
				as.integer(n.aug),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
				as.integer(maxaug),
				as.integer(maxiaug),
				as.double(neglect),
				as.integer(chr)
			    )
		n.ind = result[[7]]
		n.aug = result[[8]]	
		markONchr <- 0
		markdone <- 0
		for(c in 1:n.chr){
			#print(paste("Cromosome",c,"\n",sep=""))
			matri <- NULL
			matri2 <- NULL
			markONchr <- dim(cross$geno[[c]]$data)[2]
			#print(paste("# markers",markONchr,"\n",sep=""))
			for(j in markdone:(markdone+markONchr-1)){
			    #print(paste("Start",markdone,":End",(markdone+markONchr-1),"\n",sep=""))
				ind2 <- NULL
				pheno <- NULL
				ind2 = result[[4]][(1+(j*maxaug)):(n.aug+(j*maxaug))]
				matri <- rbind(matri,ind2)
			}
			pheno <- as.matrix(result[[5]][1:n.aug])
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
			cross$extra$Nind <- n.ind
			cross$extra$Naug <- n.aug
			cross$extra$augIND <- result[[6]][1:n.aug]
			markdone <- (markdone+markONchr)  
		}
		cross
	}else{
		stop("Error: Currently only f2 crosses can be analyzed by MQM.")
		return 
	}			
}

# end of MQMaugment.R

	
cross_good <- MQMaugment(cross)
listeria_good <- MQMaugment(listeria,maxaug=10000,maxiaug=1000,neglect=1000)