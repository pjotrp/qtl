\name{scanMQMall}
\alias{scanMQMall}

\title{ 
	scanMQMall - Multitrait analysis of a cross object 
}

\description{
	Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object.
}
\usage{
	scanMQMall(cross = NULL, cofactors = NULL, doLOG = 0,REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.clusters = 10)
}

\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{doLOG}{ Indicate if the algorithm can LOG transform traits with variation > 10.000 (DEFAULT 0 = off | 1 = on) }  
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{n.clusters}{ Number of child processes to split the job into }
}

\details{
  Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
  Returns a MQMmulti object. this object is a list of scanone objects that can be plotted using plot.scanone(result[[trait]])
}

\references{ Danny Arends <Danny.Arends@gmail.com> ; R.C. Janssen <r.c.jansen@rug.nl> ; Karl W. Brohamn <kbroman@biostat.wisc.edu>  }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
	\code{\link{multitrait}} - Dataset used in the example analysis
	\code{\link{MQMaugment}} - Augmentationroutine to estimate missing data
	\code{\link{scanMQM}} - Function called to do single trait analysis
}

\examples{
	#Doing a multitrait analysis
	library(MQMpackage)
	data(multitrait)
	result <- scanMQMall(multitrait,doLOG=1)
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
