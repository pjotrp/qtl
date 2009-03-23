#####################################################################
#
# PipelineMolgenis.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: PipelineMolgenis
#
######################################################################

######################################################################
#
# PipelineMolgenis:
#
######################################################################

PipelineMolgenis <- function(DBmarkerID=298,DBtraitID=181,saveAS="MQMResults",DBpath="http://celtic.service.rug.nl:8080/molgenis4rsandbox",Fupdate=0,...){
	cross <- CrossFromMolgenis(DBmarkerID,DBtraitID,DBpath)
	result <- scanMQMall(cross,...,plot=FALSE)
	ResultsToMolgenis(result, saveAS, DBpath, Fupdate)
}