#####################################################################
#
# scanMQMall.R
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
# scanMQMall:
#
######################################################################
setwd("D:/")

library(MQMpackage)
dyn.load("scanMQM.dll")
cross <- read.cross("csvr","","Test_Y.txt",genotype=c(1,2))
class(cross)[1] <- "f2"

scanMQMall <- function(cross= NULL,cofactors = NULL,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.clusters=10){

	
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2"){
		cat("INFO: Received an F2 cross.\n")
	}
	n.pheno <- nphe(cross)
	data <- NULL
	for(i in 1:n.pheno) {
		data[[i]] <- MQMaugment(cross,i)
	}
	if("snow" %in% installed.packages()[1:dim(installed.packages())[1]]){
		library(snow)
		cl <- makeCluster(n.clusters)
		clusterEvalQ(cl, library(MQMpackage))
		res <- parLapply(cl,data,scanMQM)
		stopCluster(cl)
	}else{
		res <- Lapply(cl,data,scanMQM)
	}
	colors <- rainbow(n.pheno)
	for(i in 1:n.pheno) {
		if(i !=1 ){
			plot(res[[i]],add=TRUE,col=colors[i])
		}else{
			plot(res[[i]],col=colors[i])
		}
	}
	res
}

# end of scanMQM.R

res <- scanMQMall(cross)
