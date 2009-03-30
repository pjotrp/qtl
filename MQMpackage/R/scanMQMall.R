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

scanMQMall <- function(cross= NULL,cofactors = NULL,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.clusters=2,doLOG=0,est.map=0,dominance=0,forceRIL=0,FF=0,plot=TRUE,verbose=TRUE){

	
	if(is.null(cross)){
		ourstop("No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		n.pheno <- nphe(cross)
		data <- NULL
	
		for(i in 1:n.pheno) {
			#we use the augmentation routinge to make a list where data[i] has trait[i] and a crossfile
			data[[i]] <- MQMaugment(cross,i,verbose=verbose)
		}
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
			res <- parLapply(cl,data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,est.map=est.map,forceRIL=forceRIL,plot=FALSE,verbose=FALSE)
			stopCluster(cl)
		}else{
			#Apply scanMQM to the data with the specified settings		
			cat("INFO: Library snow not found, so going into singlemode.\n")
			res <- lapply(data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG,est.map=est.map,forceRIL=forceRIL,plot=plot,verbose=verbose)
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
		res
	}else{
		stop("ERROR: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}
}

# end of scanMQMall.R
