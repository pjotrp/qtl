\name{scanMQMall}
\alias{scanMQMall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ scanMQMall - Multitrait analysis of a cross object }
\description{
	Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object.
}
\usage{
scanMQMall(cross = NULL, cofactors = NULL, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.clusters = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{n.clusters}{ Number of child processes to spil the job into }
}
\details{
  Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calcuation.
}
\value{
  Returns a list of scanone objects that can be plotted using plot.scanone(result[[trait]])
}
\references{ R.C. Janssen, K. Brohman, D. Arends  }
\author{ Danny Arends }
\note{ 

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	library(MQMpackage)
	data(multitrait)
	result <- scanMQMall(multitrait)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
