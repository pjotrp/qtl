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
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
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
f2perm <- MQMpermObject(f2result)			#Create a permutation object
summary(f2perm)										#What LOD score is considered significant ?
summary(f2result[[1]],alpha=0.05,perms=f2perm)		#Get significant QTL's (The first run is the original phenotype original)
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
