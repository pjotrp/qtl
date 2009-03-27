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

PipelineMolgenis <- function(DBmarkerID,DBtraitID,name="MQMResults",DBpath,...){
	cat("------------------------------------------------------------------\n")
	cat("Starting Molgenis <-> MQM <-> Molgenis automated pipeline\n")
	cat("Please be patient while all data is retrieved.\n")
	cat("INFO:Server:",DBpath,"\n")
	cat("INFO:Genotype info-tableID:",DBmarkerID," (DBmarkerID)\n")
	cat("INFO:Phenotype values-tableID:",DBtraitID," (DBtraitID)\n")
	cat("INFO:Results will be stored in a table named:",name,"\n")	
	cat("------------------------------------------------------------------\n")
	all_data <- CrossFromMolgenis(DBmarkerID,DBtraitID,trait=0,DBpath)
	all_data <- fill.geno(all_data)
	num_traits <- nphe(all_data)
	cat("------------------------------------------------------------------\n")
	cat("INFO data retrieval finished\n")
	cat("------------------------------------------------------------------\n")
	for(x in 1:num_traits){
		cat("\n\n------------------------------------------------------------------\n")
		cat("Trait",x,"/",num_traits,":",names(all_data$pheno)[1],"\n")
		cat("------------------------------------------------------------------\n")
		cat("INFO:Starting Stage1: scanning for QTL's\n")
		cat("------------------------------------------------------------------\n")
			result <- scanMQM(all_data,pheno.col=x,plot=T)
		cat("------------------------------------------------------------------\n")			
		cat("INFO:Finished with Stage1\n")	
		cat("INFO:Stage2: Storing calculated QTL's to Molgenis\n")
		cat("------------------------------------------------------------------\n")
			ResultsToMolgenis(result, name,(x-1),DBpath, Fupdate)
		cat("------------------------------------------------------------------\n")			
		cat("INFO:Finished with Stage2\n")
		cat("------------------------------------------------------------------\n")	
	}
}