\name{plot.MQMboot}
\alias{plot.MQMboot}

\title{ plot.MQMall - Plotting routine to display the results from a bootstrapMQM scan  }

\description{
	Plotting routine to display the results from a bootstrapMQM scan 
}
\usage{
plot.MQMboot(result = NULL)
}

\arguments{
  \item{result}{ Result object from scanMQMall }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
	No return, plotting routine
}
\author{ Rutger Brouwer \email{R.W.W.Brouwer@rug.nl}Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 	
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
#Simulated F2 Population
	library(MQMpackage)
	f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
	data(map10)                                             # Mouse genome
	f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
	f2result <- bootstrapMQM(f2cross)                       # Do a MQM scan of the genome
	plot.MQMboot(f2result)
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
