#####################################################################
#
# bootstrapMQM.R
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
# bootstrapMQM:
#
######################################################################

#setwd("D:/")
#library(qtl)
#library(MQMpackage)
#cross <- read.cross("csv","","Test.csv")

bootstrapMQM <- function(cross= NULL,cofactors = NULL,Phenot=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.run=10,file="MQM_output.txt",doLOG=0,reestimate=0,dominance=0,n.clusters=2,parametric=0){
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		if(!parametric){
			cat("Shuffleling traits between individuals.\n")
		}else{
			cat("Parametric bootstrapping\nCalculating new traits for each individual.\n")
		}
		n.pheno <- nphe(cross)
		data <- NULL
		#Set the Phenotype under intrest as the first
		cross$pheno[[1]] <- cross$pheno[[Phenot]]
		data[[1]] <- cross
		for(i in 2:n.run) {
			#we use the augmentation routinge to make a list where data[i] has trait[i] and a crossfile
			if(!parametric){
				neworder <- sample(nind(cross))			
				cross$pheno[[1]] <- cross$pheno[[1]][neworder]
			}else{
				
				pheno <- cross$pheno[[1]]
				variance <- var(pheno,na.rm = TRUE)
				for(j in 1:nind(cross)) {
					cross$pheno[[1]][j] <- variance^0.5*runif(1)
				}
			}
			data[[i]] <- cross
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
		if("snow" %in% installed.packages()[1:dim(installed.packages())[1]]){
			cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
			library(snow)
			cl <- makeCluster(n.clusters)
			clusterEvalQ(cl, library(MQMpackage))
			res <- parLapply(cl,data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,reestimate=reestimate)
			stopCluster(cl)
		}else{
			cat("INFO: Library snow not found, so going into singlemode.\n")
			res <- lapply(data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,reestimate=reestimate)
		}
		colors <- rainbow(n.pheno)
		for(i in 1:n.pheno) {
			if(i !=1 ){
				plot(res[[i]],add=TRUE,col=colors[i])
			}else{
				plot(res[[i]],col=colors[i])
			}
		}
		class(res) <- c(class(res),"MQMmulti")
		res
	}else{
		stop("Error: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

#result <- bootstrapMQM(cross)

# end of bootstrapMQM.R
