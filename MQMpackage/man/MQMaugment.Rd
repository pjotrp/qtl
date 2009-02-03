\name{MQMaugment}
\alias{MQMaugment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MQMaugment - Missing data augmentation routine for MQM }
\description{
  
}
\usage{
MQMaugment(cross = NULL, maxaug = 10000, maxiaug = 100, neglect = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{maxaug}{ Maximum size after augmentation }
  \item{maxiaug}{ Maximum augmentations per individual }
  \item{neglect}{ Possible genotypes to considered have a higher prob then 1:Neglect }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Cross object from R/QTL
}
\references{ ~put references to the literature/web site here ~ }
\author{ Danny Arends }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	library(MQMpackage)
	qtl <- c(3,15,3,7)												# QTL at chromosome 3
	data(map10)														# Mouse genome
	cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
	cross_good <- MQMaugment(cross)									# Augmentation of the data
	result <- scanMQM(cross_good)									# Do a MQM scan of the genome
	plot(result)													# Plot the results of the genome scan
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
