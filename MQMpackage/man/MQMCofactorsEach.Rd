\name{MQMCofactorsEach}
\alias{MQMCofactorsEach}

\title{ 
	MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM 
}
\description{
  Return a cofactorslist formatted to be used with the scanMQM function
}
\usage{
MQMCofactorsEach(cross = NULL, each = 3)
}
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{each}{ Every 'each' marker will be used as a cofactor }
}
\details{
	Not available yet
}
\value{
	An formatted list of cofactors to be used with the scanMQM algorithm
}
\references{ Danny Arends \email{Danny.Arends@gmail.com} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
}
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- MQMCofactorsEach(cross,3)					# Set Cofactors each third marker
result <- scanMQM(cross,a)					# Do a backward model selection
plot(result)							# Plot the results of the genome scan (with the best model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
