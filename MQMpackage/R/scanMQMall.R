#####################################################################
#
# scanMQMall.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: scanMQM
#
######################################################################

######################################################################
#
# scanMQMall: Contains scanMQMall routine and the plot.MQMall routine
#
######################################################################

scanMQMall <- function(cross= NULL,cofactors = NULL,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.clusters=2,FF=0,plot=TRUE,verbose=TRUE,...){

	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		cat("------------------------------------------------------------------\n")
		cat("Starting MQM multitrait analysis\n")
		cat("------------------------------------------------------------------\n")
		start <- proc.time()
		n.pheno <- nphe(cross)
		SUM <- NULL
		AVG <- NULL
		all_data <- fill.geno(cross)
		#Some tests from scanMQM repeated here so they are not hidden when using snow
		if((step.min+step.size) > step.max){
				ourstop("Surrent Step settings (step.min/step.max) would crash the algorithm")
		}
		if(step.min>0){
				ourstop("Step.min needs to be smaller than 0")
		}		
		if(step.size < 1){
				ourstop("Step.size needs to be larger than 1")
		}
		#TEST FOR SNOW CAPABILITIES
		if("snow" %in% installed.packages()[1:dim(installed.packages())[1]]){
			cat("INFO: Library snow found using ",n.clusters," Cores/CPU's/PC's for calculation.\n")
			library(snow)
			cl <- makeCluster(n.clusters)
				clusterEvalQ(cl, library(MQMpackage))
				res <- parLapply(cl,1:n.pheno, snowCoreALL,all_data=all_data,cofactors=cofactors,...)
			stopCluster(cl)
		}else{
			cat("INFO: Library snow not found, so going into singlemode.\n")
			res <- lapply(1:n.pheno, snowCoreALL,all_data=all_data,cofactors=cofactors,...)
		}
		if(FF){
			cat(rownames(res[[1]]),"\n",res[[1]][,1],"\n",res[[1]][,2],"\n",file="out.frank")
			for(i in 1:length(res)){
				cat("INFO: Saving trait",i,"in frankformat\n")
				qtl <- res[[i]]
				cat(colnames(qtl)[3],qtl[,3],"\n",file="out.frank",append = T)
			}
		}
		#Return the results
		class(res) <- c(class(res),"MQMmulti")
		#All done now plot the results
		if(plot){
			plot.MQMnice(res)
		}
		end <- proc.time()
		SUM <- (end-start)[3]
		AVG <- SUM/n.pheno	
		cat("------------------------------------------------------------------\n")
		cat("INFO: Elapsed time:",(SUM%/%3600),":",(SUM%%3600)%/%60,":",round(SUM%%60, digits=0),"(Hour:Min:Sec)\n")		
		cat("INFO: Average time per trait:",round(AVG, digits=3),"seconds\n")
		cat("------------------------------------------------------------------\n")	
		res
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

snowCoreALL <- function(x,all_data,cofactors,...){
	b <- proc.time()
	num_traits <- nphe(all_data)
	cat("------------------------------------------------------------------\n")
	cat("INFO: Starting analysis of trait (",x,"/",num_traits,")\n")
	cat("------------------------------------------------------------------\n")
	result <- scanMQM(all_data,cofactors=cofactors,pheno.col=x,plot=F,verbose=T,...)
	e <- proc.time()
	cat("------------------------------------------------------------------\n")
	cat("INFO: Done with the analysis of trait (",x,"/",num_traits,")\n")	
	cat("INFO: Calculation of trait",x,"took:",round((e-b)[3], digits=3)," seconds\n")
	cat("------------------------------------------------------------------\n")
	result
}


# end of scanMQMall.R
