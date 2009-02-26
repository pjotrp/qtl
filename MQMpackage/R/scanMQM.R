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
# scanMQM: main scanning function to the MQMpackage
#
######################################################################
#setwd("D:/")
#library(qtl)
#dyn.load("scanMQM.dll")
#cross <- read.cross("csv","","Test.csv")

# Smap <- sim.map(len=rep(200,20), n.mar=10)
# Bmap <- sim.map(len=rep(200,20), n.mar=1000)

scanMQM <- function(cross= NULL,cofactors = NULL,Phenot=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,file="MQM_output.txt",doLOG=0,reestimate=0,dominance=0){
    library(qtl)
	n.run=0
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		if(class(cross)[1] == "f2"){
			ctype = 1
		}
		if(class(cross)[1] == "bc"){
			ctype = 2
		}
		if(class(cross)[1] == "riself"){
			ctype = 3
		#	stop("Somethings still wrong in the algorithm, please analyse RIL as BC.")
		}
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
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
		if(n.run !=0){
			stop("Error: # of runs should be positive and < 10000.")
		}
		#CHECK if the phenotype exists
		if(Phenot != 1){
			cat("INFO: Selected phenotype ",Phenot,".\n")
			cat("INFO: # of phenotypes in object ",nphe(cross),".\n")
			if(nphe(cross) < Phenot || Phenot < 1){
				stop("Error: No such phenotype")
			}			
		}
		pheno <- cross$pheno[[Phenot]]
		if(var(pheno,na.rm = TRUE)> 1000){
			if(doLOG == 0){
				cat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"variation:",var(pheno,na.rm = TRUE),".\n")
				warning("INFO: Perhaps we should LOG-transform this phenotype, please set parameter: doLOG=1 to correct this error")
			}
		}
		if(doLOG != 0){
				#transform the cross file
				cross <- MQMlogPheno(cross,Phenot)
				pheno <- cross$pheno[[Phenot]]
		}
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
		for(i in 1:length(pheno)) {
			if(is.na(pheno[i])){
			  cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			  dropped <- c(dropped,i) 
			  n.ind = n.ind-1
			}
		}
		#throw em out
		if(!is.null(dropped)){
			geno <- geno[-dropped,]  
			pheno <- pheno[-dropped]
		}
		
		#CHECK for previously augmented dataset
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
		
		#CHECK if we have cofactors, so we can do backward elimination
		backward <- 0;
		if(is.null(cofactors)){
			cat("INFO: No cofactors, setting cofactors to 0\n")
			cofactors = rep(0,n.mark)
		}else{
			if(length(cofactors) != n.mark){
				cat("Error: # Cofactors != # Markers\n")		
			}else{
				cat("INFO:#",length(cofactors),"Cofactors received")
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

		if((step.min+step.size) > step.max){
				stop("Error: current Step setting would crash the algorithm")
		}
		if(step.min>0){
				stop("Error: step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				stop("Error: Step.size needs to be larger than 1")
		}
		qtlAchromo <- length(seq(step.min,step.max,step.size))
		cat("INFO: Number of locations per chromosome: ",qtlAchromo, "\n")
		result <- .C("R_scanMQM",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
                as.integer(geno),
				as.integer(chr),
				as.double(dist),
				as.double(pheno),
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
				QTL=as.double(rep(0,n.chr*qtlAchromo)),
				as.integer(reestimate),
				as.integer(ctype),
				as.integer(dominance)
			    )
		# initialize output object
		qtl <- NULL
		names <- NULL
		for(i in 1:(n.chr*qtlAchromo)) {
			#Store the result in the qtl object
			qtl <- rbind(qtl,c(ceiling(i/qtlAchromo),rep(seq(step.min,step.max,step.size),n.chr)[i],result$QTL[i]))
			#make names in the form: C L
			names <- c(names,paste("C",ceiling(i/qtlAchromo),"L",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
		}
		
		rownames(qtl) <- names
		colnames(qtl) = c("chr","pos (Cm)",paste("QTL",colnames(cross$pheno)[Phenot]))
		
		#So we can use carls plotting routines
		class(qtl) <- c(class(qtl),"scanone") 
		
		cat("Saving output to file: ",file, "\n")
		write.table(qtl,file)
		#Return the results
		qtl
	
	}else{
		stop("Error: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

#res <- scanMQM(cross)
#plot(res)

#bcqtl <- c(3,15,2)                                      # QTL at chromosome 3
#data(map10)                                             # Mouse genome
#bccross <- sim.cross(map10,bcqtl,n=100,type="bc")       # Simulate a BC Cross
#bcresult <- scanMQM(bccross)                            # Do a MQM scan of the genome
#plot(bcresult)                                          # Plot the results of the genome scan

# end of scanMQM.R
