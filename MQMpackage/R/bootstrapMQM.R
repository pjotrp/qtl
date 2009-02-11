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

bootstrapMQM <- function(cross= NULL,cofactors = NULL,Phenot=1,REMLorML=0,
                    alfa=0.02,em.iter=1000,windowsize=25.0,step.size=5.0,
					step.min=-20.0,step.max=220.0,n.run=100,file="MQM_output.txt",doLOG=0){
    library(qtl)
	if(is.null(cross)){
		stop("Error: No cross file. Please supply a valid cross object.") 
	}
	if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "riself"){
		cat("INFO: Received a valid cross file type:",class(cross)[1],".\n")
		#Here we are gonna implement bootstrapping of the MQM algorithm to get estimates for a threshold	
	}else{
		stop("Error: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
	}			
}

# end of bootstrapMQM.R
