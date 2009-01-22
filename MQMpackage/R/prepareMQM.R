#####################################################################
#
# prepareMQM.R
#
# copyright (c) 2009, Danny Arends
# last modified Fep, 2009
# first written Feb, 2009
# 
# Part of the R/qtl package
# Contains: prepareMQM, printMQMin
#
######################################################################

prepareMQM <- function(cross, name){
 f1mar <- NULL
 f2mar <- t(pull.geno(cross))
 for(i in 1:dim(f2mar)[1]) {
   f1mar <- rbind(f1mar,12)
   for(j in 1:dim(f2mar)[2]) {
	if(as.character(f2mar[i,j]) == '1'){
        f2mar[i,j] <- 'A'
    }
    if(as.character(f2mar[i,j]) == '2'){
        f2mar[i,j] <- 'B'
    }
    if(as.character(f2mar[i,j]) == '3'){
        f2mar[i,j] <- 'H'
    }
   }
 }
 rownames(f1mar) = rownames(f2mar)
 print("Going to write files\n")
 filename <- paste(name,"_F2.MAR.TXT", sep="")
 write.table(f2mar, file = filename, col.names=FALSE, quote=FALSE)
 print("f2 markers writen\n")
 
 filename <- paste(name,"_F2.QUA.TXT",sep="")
 f2qua <- pull.pheno(cross)
 write.table(t(f2qua), file = filename, row.names=FALSE,col.names=FALSE, quote=FALSE)
 print("f2 genotypes writen\n")
 MQM_in <- printMQMin(cross, name)
 
}


printMQMin <- function(cross, name){
	# Printing header
	n.individuals <- nind(cross)
	n.fam <- 1
	n.mark <- sum(nmar(cross))
	info <- NULL
	info <- rbind(info,c("Nind=",n.individuals))
	info <- rbind(info,c("Nfam=",n.fam))
	info <- rbind(info,c("Nmark=",n.mark))
	
	write.table(info, file = paste(name,"_mqm_in.txt",sep=""), append = FALSE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	# Printing the markers per chromosome
	n.chr <- nchr(cross)
    map <- pull.map(cross)
	result <- NULL
	for(i in 1:n.chr) {
		n.markers <- length(map[[i]])
		for(j in 1:n.markers){
			result <- rbind(result,c(i,row.names(as.matrix(map[[i]]))[j],as.matrix(map[[i]])[j]))
		}
	}
	result
	write.table(result, file = paste(name,"_mqm_in.txt",sep=""), append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	settings <- NULL;
	settings <- rbind(settings,c("Cross=","f2"))
	settings <- rbind(settings,c("FileF1=","-"))
	settings <- rbind(settings,c("FileF2=",paste(name,"_F2.MAR.TXT",sep="")))
	settings <- rbind(settings,c("FileY=",paste(name,"_F2.QUA.TXT",sep="")))
	settings <- rbind(settings,c("Dominance=","n"))
	settings <- rbind(settings,c("RemLorML=","0"))
	settings <- rbind(settings,c("defset=","y"))
	settings <- rbind(settings,c("real_simu=","0"))
	settings <- rbind(settings,c("perm_simu=","1 0"))
	write.table(as.matrix(settings), file = paste(name,"_mqm_in.txt",sep=""), append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)	
}

GenerateTestSets <- function(){
   library(qtl)
   data(map10)
   cross <- sim.cross(map10,n=100)
   prepareMQM(cross,"T1")

   data(map10)
   cross <- sim.cross(map10,c(1,30,1,0),n=100)
   prepareMQM(cross,"T2")

   data(map10)
   cross <- sim.cross(map10,c(8,-12,1,2),n=100)
   prepareMQM(cross,"T3")

   data(map10)
   cross <- sim.cross(map10,rbind(c(1,2,1,2),c(4,23,-1,2)),n=100)
   prepareMQM(cross,"T4") 
}

readMQMout <- function(cross = NULL, file = "mqm_out.txt", plot = FALSE){
   data <-read.table(file, quote=":")
   data <-data[-dim(data)[1],]
   if(plot){
     plot(rownames(data),data[,3],xlab="markers",ylab="QTL",main="MQM", type='n')
     lines(rownames(data),data[,3],col="red",cex=0.5,pch=20)
     lines(rownames(data),data[,4],col="blue",cex=0.5,pch=20)   
   }
   data
   #should be pushed to the cross object
}