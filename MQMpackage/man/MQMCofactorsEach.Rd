\name{MQMCofactorsEach}
\alias{MQMCofactorsEach}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM }
\description{
  Return a cofactorslist formatted to be used with the scanMQM function
}
\usage{
MQMCofactorsEach(cross = NULL, each = 3)
}
%- maybe also 'usage' for other objects documented here.
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
\references{ ~put references to the literature/web site here ~ }
\author{ Danny Arends }
\note{ 
	If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- MQMCofactorsEach(cross,3)						# Set Cofactors each third marker
result <- scanMQM(cross,a)					# Do a backward model selection
plot(result)							# Plot the results of the genome scan (with the best model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
