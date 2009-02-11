\name{scanMQM}
\alias{scanMQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ scanMQM - Multiple qtl mapping for inbred crosses }
\description{
Main multiQTL scanning function for the R/qtl package (for use with F2 / BC / RIL populations)
}
\usage{
scanMQM(cross = NULL, cofactors = NULL, doLOG = 0, Phenot = 1, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.run = 0, file = "MQM_output.txt")
}
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{Phenot}{ Phenotype to do thepull from the cross object }    
  \item{doLOG}{ Indicate if the algorithm can LOG transform traits with variation > 10.000 (DEFAULT 0 = off | 1 = on) }    
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{n.run}{ When setting n.run >= 1, a simple permuation of the data will occur, giving insight into the significance of the QTL }
  \item{file}{ File to which the raw results are saved (MQM_output.txt) }
}
\details{
Using ML or Reduced ML the algorithm employs a forward- or backward selection strategy to identify QTL hotspot.
The algoreithm passes through the following stages:
Calculating relative marker positions
Calculating recombinant frequencies
IF(Missing data) THEN Fix for dataaugmentation
Rmixture function
Loglikelyhood estimated of full model
IF(Cofactors) THEN Backward elimination & mapQTL
IF(No cofactors) THEN mapQTL function
IF(Nrun) THEN Bootstrapping either mapQTL or Backward elimination

}
\value{
Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the plot.scanone function of R/QTL.
The phenotype under consideration is specified by the Phenot argument of the function.
}
\references{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{  

}
\examples{
#Simulated F2 Population
library(MQMpackage)
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- scanMQM(f2cross)				# Do a MQM scan of the genome
plot(f2result)						# Plot the results of the genome scan

#Simulated BC Population
library(MQMpackage)
bcqtl <- c(3,15,2)					# QTL at chromosome 3
data(map10)						# Mouse genome
bccross <- sim.cross(map10,bcqtl,n=100,type="bc")	# Simulate a BC Cross
bcresult <- scanMQM(bccross)				# Do a MQM scan of the genome
plot(bcresult)						# Plot the results of the genome scan

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
