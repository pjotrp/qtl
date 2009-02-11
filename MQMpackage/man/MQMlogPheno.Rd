\name{MQMlogPheno}
\alias{MQMlogPheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MQMlogPheno - Helper function to logtransform phenotypesin R/QTL cross objects }
\description{
  Helper function to logtransform phenotypes in an R/QTL cross objects cross&pheno <- log(cross$pheno)
}
\usage{
MQMlogPheno(cross = NULL, Phenot = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{Phenot}{ Phenot to logtransform (if none is specified all phenotypes are transformed) }
}
\details{
	Not available yet
}
\value{
	Cross object with log transformed phenotypes (specified bij Phenot, or all are transformed when Phenot=NULL) 
}
\references{ ~put references to the literature/web site here ~ }
\author{  Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 

}
\seealso{ 
	\code{\link{scanMQM}} - Function called to do single trait analysis
	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
\examples{
	data(multitrait)
	plot(scanMQM(MQMaugment(multitrait)))
	plot(scanMQM(MQMaugment(multitrait),doLOG=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
