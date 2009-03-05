\name{oat}
\alias{oat}

\docType{data}

\title{ 
	Example F2 Cross object from R/QTL with a single trait.
}

\description{
	Cross object from R/QTL, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details.
}
\usage{data(multitrait)}
\format{
	Cross object from R/QTL
}
\details{
	An F2 inbred line of Kanota Oat with Ogle oat
}
\source{
  1 Trait of the oat breading experiment with Kanota (K) and Ogle (O) with 50 individuals scored (with errors at) 202 markers.
  Dataset obtained from GBIC - Groningen BioInformatics Centre
}
\references{
  \item Wight et al, (2003) \emph{Genome} \bold{46},28-47
  \item O'Donoughue et al, (1995) \emph{Genome} \bold{38},368-380
}
\examples{
	library(MQMpackage)			#Load the package
	data(oat)		#Load dataset
	result <- scanMQM(oat)	#Do a basic MQMscan of the genome
}
\keyword{datasets}
