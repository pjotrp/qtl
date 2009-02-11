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
# scanMQMall: Contains scanMQMall routine and the plot.MQMall routine
#
######################################################################
scanMQMall <- function(cross= NULL,cofactors = NULL,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.clusters=2,doLOG=0){

	
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "ril"){
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		n.pheno <- nphe(cross)
		data <- NULL
	
		for(i in 1:n.pheno) {
			#we use the augmentation routinge to make a list where data[i] has trait[i] and a crossfile
			data[[i]] <- MQMaugment(cross,i)
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
			res <- parLapply(cl,data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG)
			stopCluster(cl)
		}else{
			cat("INFO: Library snow not found, so going into singlemode.\n")
			res <- lapply(data,scanMQM,step.min=step.min,step.max=step.max,alfa=alfa,em.iter=em.iter,windowsize=windowsize,REMLorML=REMLorML,cofactors=cofactors,step.size=step.size,doLOG=doLOG)
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

# end of scanMQMall

plot.MQMall <- function(cross=NULL, result = NULL, type="C",theta=30,phi=15){
	if(is.null(cross)|| is.null(result)){
		stop("Error: NO result or cross.") 
	}
	if(class(result)[2] == "MQMmulti"){
		if(type=="C"){
			c <- NULL
			for(i in 1:length(result)){
				c <- rbind(c,result[[i]][,3])
			}
			c <- t(c)
			contour(
				x=seq(1,dim(c)[1]),
				y=seq(1,dim(c)[2]),
				c,
				xlab="Markers",ylab="Trait",
				col=rainbow((max(c)/5)+25,1,1.0,0.1),
				nlevels=(max(c)/5)
			)
		}
		if(type=="I"){
			c <- NULL
			for(i in 1:length(result)){
				c <- rbind(c,result[[i]][,3])
			}
			c <- t(c)
			image(x=1:dim(c)[1],y=1:dim(c)[2],c,
				  xlab="Markers",ylab="Trait",
				  col=rainbow((max(c)/5)+25,1,1.0,0.1),
			)
		
		}
		if(type=="D"){
			c <- NULL
			for(i in 1:length(result)){
				c <- rbind(c,result[[i]][,3])
			}
			c <- t(c)
			persp(x=1:dim(c)[1],y=1:dim(c)[2],c,
				  theta = theta, phi = phi, expand = 1,
				  col="gray", xlab = "Markers", ylab = "Traits", zlab = "QTL")
		}
		if(type=="P"){
			n.pheno <- nphe(cross)
			colors <- rainbow(n.pheno)
			for(i in 1:n.pheno) {
				if(i !=1 ){
					plot(result[[i]],add=TRUE,col=colors[i])
				}else{
					plot(result[[i]],col=colors[i])
				}
			}
		}			
	}else{
		stop("Error: wrong type of result file.") 
	}
}

# end of plot.MQMall
