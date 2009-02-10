\name{MQMCofactors}
\alias{MQMCofactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MQMCofactors - Helper function to create cofactors to be used with scanMQM }
\description{
  Return a cofactorslist formatted to be used with the scanMQM function
}
\usage{
MQMCofactors(cross = NULL, cofactors = NULL, sexfactors = NULL)
}
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of integer describing which marker should be treated as a cofactor }
  \item{sexfactors}{ List of integer describing which marker should be treated as a sexfactor }
}
\details{
Not available yet
}
\value{
An formatted list of cofactors to be used with the scanMQM algorithm
}
\references{ ~put references to the literature/web site here ~ }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{  
	\code{\link{scanMQM}} - Function called to do single trait analysis
	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- MQMCofactors(cross,c(20,30,80,90),c(186,187))		# Set Cofactors on marker 20,30,80 & 90 and sexfactors on 186,187
result <- scanMQM(cross,a)					# Do a backward model selection
plot(result)							# Plot the results of the genome scan (with the best model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
