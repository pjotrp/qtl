\name{MQMpermObject}
\alias{MQMpermObject}

\title{ MQMpermObject - Function to convert MQMmulti objects into a scanoneperm object from R/QTL }
\description{
	Function to convert MQMmulti objects into a scanoneperm object from R/QTL, this allows us to use the R/QTL methods of permutation analysis
}
\usage{
MQMpermObject(MQMbootresult = NULL)
}

\arguments{
  \item{MQMbootresult}{ MQMmulti object obtained by bootstrapMQM }
}
\details{

}
\value{
	a scanoneperm object
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} }
\note{ 

}
\seealso{ 
\itemize{
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
library(MQMpackage)
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- bootstrapMQM(f2cross)			# Do a MQM scan of the genome
f2.perm <- MQMpermObject(f2result)		#Create a permutation object
summary(f2result,alpha=0.05,perms=f2.perm)	#get significant QTL's
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
