#####################################################################
#
# MQMLicence.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: MQMLicence
#
######################################################################

MQMLicence <- function(){
    text <- paste("This software is distributed under the terms of \n")
	text <- paste("Developed by: RC Janssen\nImplementation into R: Danny Arends\n")
	MQMwindow(text)
}

MQMwindow <- function(text){
    outFile <- tempfile()
    outConn <- file(outFile, open = "w")
    writeLines(paste(text, sep = ""), outConn)
    close(outConn)
    file.show(outFile, delete.file = TRUE)
}