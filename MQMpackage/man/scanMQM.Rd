\name{scanMQM}
\alias{scanMQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ scanMQM - Multiple qtl mapping for F2 crosses }
\description{
Main multiQTL scanning function for the R/qtl package (for use with F2 populations)
}
\usage{
scanMQM(cross = NULL, cofactors = NULL, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.run = 0)
}
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{Phenot}{ Phenotype to do thepull from the cross object }    
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Window size when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL. }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{n.run}{ When setting n.run >= 1, a simple permuation of the data will occur, giving insight into the significance of the QTL }
}
\details{
Using ML or Reduced ML the algorithm employs a forward- or backward selection strategy to identify QTL hotspot.
}
\value{
Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the plot.scanone function of R/QTL.
The phenotype under consideration is specified by the Phenot argument of the function.
}
\references{ R.C. Janssen, K. Brohman, D. Arends }
\author{ Danny Arends }
\note{
If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{  }
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- scanMQM(cross)				# Do a MQM scan of the genome
plot(result)						# Plot the results of the genome scan
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
