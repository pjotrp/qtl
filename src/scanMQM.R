#####################################################################
#
# scanMQM.R
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
# scanMQM:
#
######################################################################

scanMWM <- function(cross){
	result <- .C("R_scanMQM")
}

# end of scanMQM.R