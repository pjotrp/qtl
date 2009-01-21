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

prepareMQM <- function(cross){

 f2mar <- t(pull.geno(cross))
 for(i in 1:dim(f2mar)[1]) {
   for(j in 1:dim(f2mar)[2]) {
	if(as.character(f2mar[i,j]) == '1'){
        f2mar[i,j] <- 11
    }
    if(as.character(f2mar[i,j]) == '2'){
        f2mar[i,j] <- 22
    }
    if(as.character(f2mar[i,j]) == '3'){
        f2mar[i,j] <- 12
    }
   }
 }
 write.table(f2mar, file = "GEN_F2.MAR.TXT", col.names=FALSE, quote=FALSE)
 f2qua <- pull.pheno(cross)
 write.table(f2qua, file = "GEN_F2.QUA.TXT", row.names=FALSE,col.names=FALSE, quote=FALSE)
 MQM_in <- printMQMin(cross)
 
}


printMQMin <- function(cross){
	# Printing header
	n.individuals <- nind(cross)
	n.fam <- 1
	n.mark <- sum(nmar(cross))
	info <- NULL
	info <- rbind(info,c("Nind=",n.individuals))
	info <- rbind(info,c("Nfam=",n.fam))
	info <- rbind(info,c("Nmark=",n.mark))
	
	write.table(info, file = "GEN_MQM_IN.TXT", append = FALSE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
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
	write.table(result, file = "GEN_MQM_IN.TXT", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)
	
	settings <- NULL;
	settings <- rbind(settings,c("Cross=","f2"))
	settings <- rbind(settings,c("FileF1=","XXX"))
	settings <- rbind(settings,c("FileF2=","GEN_F2.MAR.TXT"))
	settings <- rbind(settings,c("FileY=","GEN_F2.QUA.TXT"))
	settings <- rbind(settings,c("Dominance=","0"))
	settings <- rbind(settings,c("RemLorML=","0"))
	settings <- rbind(settings,c("defset=","1"))
	settings <- rbind(settings,c("real_simu=","0"))
	settings <- rbind(settings,c("perm_simu=","0"))
	write.table(as.matrix(settings), file = "GEN_MQM_IN.TXT", append = TRUE, sep = " ",col.names=FALSE,row.names=FALSE, quote=FALSE)	
}