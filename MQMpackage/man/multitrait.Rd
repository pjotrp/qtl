\name{multitrait}
\alias{multitrait}

\docType{data}

\title{ 
	Example Cross object from R/QTL with multiple traits.
}

\description{
	Cross object from R/QTL, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details.
}
\usage{data(multitrait)}
\format{
	Cross object from R/QTL
}
\details{
	An f2 recombinant inbred line
}
\source{
  GBIC - Groningen BioInformatics Centre
}
\references{
  ~~ possibly secondary sources and usages ~~
}
\examples{
	library(MQMpackage)
	data(multitrait)
	result <- scanMQMall(multitrait,doLOG=1)
}
\keyword{datasets}
