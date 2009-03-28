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
	cat("INFO:Server:",DBpath,"\n")
	cat("INFO:Genotype info-tableID:",DBmarkerID," (DBmarkerID)\n")
	cat("INFO:Phenotype values-tableID:",DBtraitID," (DBtraitID)\n")
	cat("INFO:Results will be stored in a table named:",name,"\n")	
	cat("------------------------------------------------------------------\n")
	cat("Starting data retrieval.\n")	
	cat("Please be patient while all data is retrieved.\n")	
	cat("------------------------------------------------------------------\n")
	start <- proc.time()
	all_data <- CrossFromMolgenis(DBmarkerID,DBtraitID,trait=0,DBpath,verbose=F)
	all_data <- fill.geno(all_data)
	num_traits <- nphe(all_data)
	end <- proc.time()		
	cat("------------------------------------------------------------------\n")
	cat("INFO data retrieval finished in ",round((end-start)[3], digits = 3)," seconds\n")
	cat("------------------------------------------------------------------\n")
	SUM <- 0
	AVG <- 0
	LEFT <- 0
	for(x in 1:num_traits){
		start <- proc.time()
		cat("\n\n------------------------------------------------------------------\n")
		cat("Trait",x,"/",num_traits,":",names(all_data$pheno)[x],"\n")
		cat("------------------------------------------------------------------\n")
		cat("INFO:Starting Stage1: scanning for QTL's\n")
		cat("------------------------------------------------------------------\n")
			result <- scanMQM(all_data,pheno.col=x,plot=T,verbose=F,...)
		cat("------------------------------------------------------------------\n")			
		cat("INFO:Finished with Stage1\n")	
		cat("INFO:Stage2: Storing calculated QTL's to Molgenis\n")
		cat("------------------------------------------------------------------\n")
			ResultsToMolgenis(result, name,(x-1),DBpath, Fupdate,verbose=F)
		end <- proc.time()		
		SUM <- SUM + (end-start)[3]
		AVG <- SUM/x
		LEFT <- AVG*(num_traits-x)
		cat("------------------------------------------------------------------\n")
		cat("INFO:Finished with Stage\n")
		cat("INFO:Calculation of trait ",x," took: ",round((end-start)[3], digits=3)," seconds\n")
		cat("INFO:Elapsed time:",round(SUM, digits=3),"seconds\n")
		cat("INFO:Average time per trait:",round(AVG, digits=3),"seconds\n")
		cat("INFO:estimated time left:",round(LEFT, digits=3),"seconds\n")
		cat("------------------------------------------------------------------\n")	
	}
}