\name{scanMQM}
\alias{scanMQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ scanMQM - Multiple qtl mapping }
\description{
  Function for the R/qtl package
}
\usage{
scanMQM(cross = NULL, cofactors = NULL, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors (0:no cofactor, 1 cofactor, 2 sexfactor }
  \item{REMLorML}{ Reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ P-value (0<= alfa <= 1 }
  \item{em.iter}{ Max itterations of the EM algorithm }
  \item{windowsize}{ Window size when mappingQTL locations }
  \item{step.size}{ Stepsize used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  No return value, the cross objkect will be updated with the new information
}
\references{ R.C. Janssen, C. Brohman, Danny Arends }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{  }
\examples{
	library(qtl)
	qtl <- c(3,15,+3,+7)							# QTL at chromosome 3
	map10 <- data(map10)							# Mouse genome
	cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
	result <- scanMQM(cross)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
