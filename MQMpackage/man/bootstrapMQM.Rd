\name{bootstrapMQM}
\alias{bootstrapMQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ bootstrapMQM - Bootstrapping approach to obtain estimates of QTL significance }
\description{
   Bootstrapping approach to obtain estimates of QTL significance
}
\usage{
bootstrapMQM(cross = NULL, cofactors = NULL, Phenot = 1, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.run = 10, file = "MQM_output.txt", doLOG = 0, reestimate = 0, dominance = 0, n.clusters = 2, parametric = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{n.run}{ Number of bootstrapping runs (DEFAULT=10) }
  \item{parametric}{ What kind of bootstrapping should be done ( 0= Random Shuffle (DEFAULT)/ 1= Parametric bootstrapping ) }  
  \item{Phenot}{ Phenotype to do thepull from the cross object }
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{doLOG}{ Indicate if the algorithm can LOG transform traits with variation > 10.000 (DEFAULT 0 = off | 1 = on) }
  \item{reestimate}{ Reestimation of the genetic map (est.rf in R/QTL) }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }
  \item{n.clusters}{ When SNOW is installed use this to allocate more core's }
  
}
\details{
  scanMQM to scan for QTL's in bootstrapped data. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}
\value{
	Returns a MQMmulti object. this object is a list of scanone objects that can be plotted using plot.scanone(result[[trait]])
}
\references{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{
	\itemize{
	\item \code{\link{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link{scanMQMall}} - Function called to do single trait analysis
	}
}
\examples{
#Simulated F2 Population
library(MQMpackage)
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- bootstrapMQM(f2cross)			# Do a MQM scan of the genome
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
