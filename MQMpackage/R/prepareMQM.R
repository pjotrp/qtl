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

prepareMQM <- function(cross, name,cofactors=NULL,dominance='n',RemLorML=0){
#print files needed by the MQM algorithm (Genotype & Phenotype, also calls 
#the construction function of the run-input file)

 f1mar <- NULL
 f2mar <- t(pull.geno(cross))
 n.ind <- nind(cross)
 for(i in 1:dim(f2mar)[1]) {
   f1mar <- rbind(f1mar,12)
   for(j in 1:dim(f2mar)[2]) {
    if(is.na(f2mar[i,j])){
	     f2mar[i,j] <- '-'
	}
	if(as.character(f2mar[i,j]) == '1'){
        f2mar[i,j] <- 'A'
    }
    if(as.character(f2mar[i,j]) == '2'){
        f2mar[i,j] <- 'H'
    }
    if(as.character(f2mar[i,j]) == '3'){
        f2mar[i,j] <- 'B'
    }
    if(as.character(f2mar[i,j]) == '4'){
        f2mar[i,j] <- 'C'
    }
    if(as.character(f2mar[i,j]) == '5'){
        f2mar[i,j] <- 'D'
    }	
	}
 }
 rownames(f1mar) = rownames(f2mar)
 pheno <- cross$pheno
 #check for missing phenotypes
 dropped <- NULL
 for(i in 1:dim(pheno)[1]) {
		if(is.na(pheno[i,1])){
			cat("INFO: Dropped individual ",i," with missing phenotype.\n")
			dropped <- c(dropped,i) 
			n.ind = n.ind-1
		}
	}
 #throw em out
    if(!is.null(dropped)){
			f2mar<- f2mar[,-dropped]  
			pheno <- pheno[-dropped,]
		}
 filename <- paste(name,"_F2.MAR.TXT", sep="")
 write.table(f2mar, file = filename, col.names=FALSE, quote=FALSE)

 filename <- paste(name,"_F2.QUA.TXT",sep="")
 write.table(t(pheno), file = filename, row.names=FALSE,col.names=FALSE, quote=FALSE)
 MQM_in <- printMQMin(cross,name,cofactors,dominance,RemLorML) 
}

printMQMin <- function(cross,name,cofactors=NULL,dominance='n',RemLorML=0){
#print mqm_in.txt file needed by the MQM algorithm
	pheno <- cross$pheno	
	n.individuals <- nind(cross)
	for(i in 1:dim(pheno)[1]) {
		if(is.na(pheno[i,1])){
			n.individuals = n.individuals-1
		}
	}	
	n.fam <- 1
	n.mark <- sum(nmar(cross))
	info <- NULL
	info <- rbind(info,c("Nind=",n.individuals))
	info <- rbind(info,c("Nfam=",n.fam))
	info <- rbind(info,c("Nmark=",n.mark))
	
	write.table(info, file = "mqm_in.txt", append = FALSE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	# Printing the markers per chromosome
	n.chr <- nchr(cross)
	cur.marker <- 0
    map <- pull.map(cross)
	result <- NULL
	for(i in 1:n.chr) {
		n.markers <- length(map[[i]])
		for(j in 1:n.markers){
		    if(cur.marker %in% cofactors){
				result <- rbind(result,c(i,row.names(as.matrix(map[[i]]))[j],as.matrix(map[[i]])[j],"*"))
			}else{
			    result <- rbind(result,c(i,row.names(as.matrix(map[[i]]))[j],as.matrix(map[[i]])[j]," "))
			}
			cur.marker <- cur.marker+1
		}
	}
	write.table(result, file = "mqm_in.txt", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	settings <- NULL;
	settings <- rbind(settings,c("Cross=","f2"))
	settings <- rbind(settings,c("FileF1=","-"))
	settings <- rbind(settings,c("FileF2=",paste(name,"_F2.MAR.TXT",sep="")))
	settings <- rbind(settings,c("FileY=",paste(name,"_F2.QUA.TXT",sep="")))
	settings <- rbind(settings,c("Dominance=",dominance))
	settings <- rbind(settings,c("RemLorML=",RemLorML))
	settings <- rbind(settings,c("defset=","y"))
	settings <- rbind(settings,c("real_simu=","0"))
	settings <- rbind(settings,c("perm_simu=","1 0"))
	write.table(as.matrix(settings), file = "mqm_in.txt", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)	
}

RunTest <- function(testset = "T1", exe="V_1.exe"){
#Runs an executable versus a, and compares the output to the output from runPrelimTest 
 
 shell("del mqm_out.txt", intern=TRUE) 
 #read old output
 v0.output <- read.table(paste("V_0/",testset,"_OUT",sep=""),comment.char = ":")
 #setup the In-File for testing
 outcome <- shell(paste("copy V_0\\",testset,"_mqm_in.txt mqm_in.txt",sep=""), intern=TRUE)
 outcome <- shell(exe, intern=TRUE)
 vNew_out <- readMQMout()
 error <- 0
 for(i in 1:dim(v0.output)[1]) {
   for(j in 1:dim(v0.output)[2]) {
    if(abs(v0.output[i,j] - vNew_out[i,j]) > 0.01){
	#print(paste("WARNING No match between old and new version at: ",i,j,"",sep=""))
	error <- 1;
	}
  }
 }
 if(error == 0){
   print(paste("No error detected between V_0.exe and ",exe,sep=""))
 }else{
   print(paste("Errors detected between V_0.exe and ",exe,sep=""))
 }
}

RunPrelimTest <- function(testset = "T1",qtl = NULL,cofactors=NULL,dominance='n',RemLorML=0){
#Runs MQM (original version on a testset)
#Places files in an "V_0" directory (which should be setup in advance by the user

   library(qtl)
   setwd("D:/MQM/compiled")
   library(qtl)
   shell("del mqm_out.txt", intern=TRUE)    # make sure the output from the previous runs are gone
   data(map10)
   cross <- sim.cross(map10,qtl,n=100)
   prepareMQM(cross,testset,cofactors,dominance,RemLorML)
   outcome <- shell("V_0.exe", intern=TRUE)
   write.table(readMQMout(), file = paste("V_0/",testset,"_OUT",sep=""), row.names=FALSE,col.names=FALSE, quote=FALSE)
   shell(paste("copy mqm_in.txt V_0\\",testset,"_mqm_in.txt",sep=""), intern=TRUE)   #Stores the runFILE & output to a V_0 directory
   shell("del mqm_out.txt", intern=TRUE)   											#make sure
   shell("del mqm_in.txt", intern=TRUE)    											# also replace the input file (for the new run
}


GenerateTestSets <- function(){
#Run once version to setup some testfiles
#Executes "V_0.exe" to get output for each testfile (so we can compare newer versions)

   RunPrelimTest(testset="T1")
   RunPrelimTest(testset="T1c",cofactors=c(10,50,100))   
   RunPrelimTest(testset="T2",qtl=c(1,30,1,0))
   RunPrelimTest(testset="T2c",qtl=c(1,30,1,0),cofactors=c(10,50,100))
   RunPrelimTest(testset="T3",qtl=c(19,30,0,1))
   RunPrelimTest(testset="T4",qtl=c(19,500,0,1))
   RunPrelimTest(testset="T5",qtl=rbind(c(19,5,0,1),c(19,45,0,1)))
   RunPrelimTest(testset="T6",qtl=rbind(c(19,5,0,1),c(19,45,0,-1)))

}

readMQMout <- function(cross = NULL, file = "mqm_out.txt", plot = FALSE,chr = 1){
#reads the output from the MQM algorithm
   data <-read.table(file, quote=":")
   data <-data[-dim(data)[1],]
   if(plot){
       plot(rownames(data),data[,1]/1.5,xlab="Markers",ylab="QTL",main="MQM", type='n')
       lines(rownames(data),data[,3],col="red",cex=0.5,pch=20)
       lines(rownames(data),data[,4],col="blue",cex=0.5,pch=20)  
   }
   rownames(data) <- paste(data[,1],data[,2])
   data <- data[,-4]
   colnames(data) <- c("chr","pos","lod")
   data
   #should be pushed to the cross object
}