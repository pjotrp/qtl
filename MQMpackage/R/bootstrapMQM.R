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
# bootstrapMQM: Shuffles phenotype or does parametric bootstrapping of scanMQM
#
######################################################################

#setwd("D:/")
#library(qtl)
#library(MQMpackage)
#cross <- read.cross("csv","","Test.csv")

bootstrapMQM <- function(cross= NULL,cofactors = NULL,Phenot=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.run=10,file="MQM_output.txt",doLOG=0,reestimate=0,dominance=0,n.clusters=2,parametric=0)
{
	
	if(is.null(cross)){
		stop("ERROR: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		#Echo back the cross type
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		if(!parametric){
			cat("INFO: Shuffleling traits between individuals.\n")
		}else{
			cat("INFO: Parametric bootstrapping\nINFO: Calculating new traits for each individual.\n")
		}
		n.pheno <- nphe(cross)
		data <- NULL
		#Set the Phenotype under intrest as the first
		cross$pheno[[1]] <- cross$pheno[[Phenot]]
		#Set the first run to not do anything with the data
		data[[1]] <- cross
		for(i in 2:n.run) {
			#For each run either shuffle the trait-values or calculate 
			#new random ones per individuals based on the variance
			if(!parametric){
				neworder <- sample(nind(cross))			
				cross$pheno[[1]] <- cross$pheno[[1]][neworder]
			}else{
				
				pheno <- cross$pheno[[1]]
				variance <- var(pheno,na.rm = TRUE)
				for(j in 1:nind(cross)) {
					cross$pheno[[1]][j] <- runif(1)*(variance^0.5)
				}
			}
			#Set the data for this run to be the new cross object
			data[[i]] <- cross
		}
		#Some tests from scanMQM repeated here so they are not hidden when using snow
		if((step.min+step.size) > step.max){
				stop("ERROR: current Step setting would crash the algorithm")
		}
		if(step.min>0){
				stop("ERROR: step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				stop("ERROR: Step.size needs to be larger than 1")
		}
		#TEST FOR SNOW CAPABILITIES
		if("snow" %in% installed.packages()[1:dim(installed.packages())[1]]){
			cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
			library(snow)
			cl <- makeCluster(n.clusters)
			clusterEvalQ(cl, library(MQMpackage))
			res <- parLapply(cl,data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,reestimate=reestimate,plot=FALSE)
			stopCluster(cl)
		}else{
			#Apply scanMQM to the data with the specified settings
			cat("INFO: Library snow not found, so going into singlemode.\n")
			res <- lapply(data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,reestimate=reestimate,plot=FALSE)
		}
		
		#Set the class of the result to MQMmulti (so we can use our plotting routines)
		class(res) <- c(class(res),"MQMmulti")
		#All done now plot the results
		plot.MQMall(res,"P")
		#Return the results	
		res
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

#result <- bootstrapMQM(cross)

# end of bootstrapMQM.R
