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
#library(MQMpackage)
#dyn.load("scanMQM.dll")
#cross <- read.cross("csv","","Test.csv")
#cof <- MQMCofactorsEach(cross,10)
	
	
scanMQM <- function(cross= NULL,cofactors = NULL,pheno.col=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,file="MQM_output.txt",doLOG=0,est.map=0,dominance=0,plot=TRUE){
    library(qtl)
	n.run=0
	if(is.null(cross)){
		stop("ERROR: No cross file. Please supply a valid cross object.") 
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
			stop("ERROR: Alfa must be between 0 and 1.")
		}
		if(n.run !=0){
			stop("ERROR: # of runs should be positive and < 10000.")
		}
		#CHECK if the phenotype exists
		if(pheno.col != 1){
			cat("INFO: Selected phenotype ",pheno.col,".\n")
			cat("INFO: Number of phenotypes in object ",nphe(cross),".\n")
			if(nphe(cross) < pheno.col || pheno.col < 1){
				stop("ERROR: No such phenotype in cross object.\n")
			}			
		}
		pheno <- cross$pheno[[pheno.col]]
		if(var(pheno,na.rm = TRUE)> 1000){
			if(doLOG == 0){
				cat("INFO: Before LOG transformation Mean:",mean(pheno,na.rm = TRUE),"variation:",var(pheno,na.rm = TRUE),".\n")
				warning("INFO: Perhaps we should LOG-transform this phenotype, please set parameter: doLOG=1 to correct this error")
			}
		}
		if(doLOG != 0){
				#transform the cross file
				cross <- MQMlogPheno(cross,pheno.col)
				pheno <- cross$pheno[[pheno.col]]
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
				cat("ERROR: # Cofactors != # Markers\n")		
			}else{
				cat("INFO:",length(cofactors),"Cofactors received to be analyzed\n")
				if(sum(cofactors) > 0){
					cat("INFO: Doing backward elimination of selected cofactors.\n")
					backward <- 1;
					n.run <- 0;
				}else{
					backward <- 0;
					stop("ERROR: Are u trying to give an empty cofactor list ???")
				}
			}
		}

		if((step.min+step.size) > step.max){
				stop("ERROR: Current Step setting would crash the algorithm")
		}
		if(step.min>0){
				stop("ERROR: step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				stop("ERROR: Step.size needs to be larger than 1")
		}
		qtlAchromo <- length(seq(step.min,step.max,step.size))
		cat("INFO: Number of locations per chromosome: ",qtlAchromo, "\n")
		result <- .C("R_scanMQM",
				as.integer(n.ind),
                as.integer(n.mark),
				as.integer(1),    # 1 phenotype
                as.integer(geno),
				as.integer(chr),
				DIST=as.double(dist),
				as.double(pheno),
				COF=as.integer(cofactors),
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
				QTL=as.double(rep(0,2*n.chr*qtlAchromo)),
				as.integer(est.map),
				as.integer(ctype),
				as.integer(dominance)
			    )
		# initialize output object
		qtl <- NULL
		info <- NULL
		names <- NULL
		for(i in 1:(n.chr*qtlAchromo)) {
			#Store the result in the qtl object
			qtl <- rbind(qtl,c(ceiling(i/qtlAchromo),rep(seq(step.min,step.max,step.size),n.chr)[i],result$QTL[i]))
			info <- rbind(info,result$QTL[(n.chr*qtlAchromo)+i])
			#make names in the form: C L
			names <- c(names,paste("C",ceiling(i/qtlAchromo),"L",rep(seq(step.min,step.max,step.size),n.chr)[i],sep=""))
		}
		if(plot){
			if(est.map && backward){
				op <- par(mfrow = c(3,1))
			}else{
				if(est.map || backward){
					op <- par(mfrow = c(2,1))
				}else{
					op <- par(mfrow = c(1,1))
				}
			}
			if(est.map){
				new_map <- pull.map(cross)
				aa <- nmar(cross)
				sum <- 1
				for(i in 1:length(aa)) {
					for(j in 1:aa[[i]]) {
						new_map[[i]][j] <- result$DIST[sum]
						sum <- sum+1
					}
				}
				cat("INFO: Viewing the user supplied map versus genetic map used during analysis.\n")
				plot.map(pull.map(cross), new_map,main="Supplied map versus re-estimated map")
			}
			if(backward){
				if(!est.map){
					new_map <- pull.map(cross)
				}
				aa <- nmar(cross)			
				sum <- 1
				qc <- NULL
				qp <- NULL
				qn <- NULL
				for(i in 1:length(aa)) {
					for(j in 1:aa[[i]]) {
						#cat("INFO ",sum," ResultCOF:",result$COF[sum],"\n")
						if(result$COF[sum] != 48){
							cat("MODEL: Marker",sum,"from model found, CHR=",i,",POSITION=",as.double(unlist(new_map)[sum])," Cm\n")
							qc <- c(qc, as.character(i))
							qp <- c(qp, as.double(unlist(new_map)[sum]))
							qn <- c(qn, substr(names(unlist(new_map))[sum],3,nchar(names(unlist(new_map))[sum])))
						}
						sum <- sum+1
					}
				}
				why <- sim.geno(cross)
				if(!is.null(qc)){
					qtlplot <- makeqtl(why, qc, qp, qn, what="draws")
					plot(qtlplot)
				}
			}
		}
		rownames(qtl) <- names
		qtl <- cbind(qtl,1/(min(info))*(info-min(info)))
		qtl <- cbind(qtl,1/(min(info))*(info-min(info))*qtl[,3])
		colnames(qtl) = c("chr","pos (Cm)",paste("QTL",colnames(cross$pheno)[pheno.col]),"Info","QTL*INFO")
		#So we can use carls plotting routines
		qtl <- as.data.frame(qtl)
		class(qtl) <- c("scanone",class(qtl)) 
		
		cat("INFO: Saving output to file: ",file, "\n")
		write.table(qtl,file)
		#Reset plotting and return the results
		if(plot){
			info_c <- qtl
			#Check for error in the information content
			e <- NULL
			for(i in 1:ncol(qtl)){
				if(is.na(info_c[i,5])){
					e<- 1
				}
				if(is.infinite(info_c[i,5])){
					e<- 1
				}
				if(is.null(info_c[i,5])){
					e<- 1
				}
			}
			#No error plot 2
			if(!e){
				info_c[,3]<- info_c[,5]
				plot(qtl,info_c,lwd=1)
				labels <- c(paste("QTL",colnames(cross$pheno)[pheno.col]),"QTL * Info")
				legend("topright", labels,col=c("black","blue"),lty=c(1,1))
			}else{
				plot(qtl,lwd=1)
			}
		}
		#Reset the plotting window to contain 1 plot
		op <- par(mfrow = c(1,1))
		qtl
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

#res <- scanMQM(cross,cof)
#plot(res)


# end of scanMQM.R
