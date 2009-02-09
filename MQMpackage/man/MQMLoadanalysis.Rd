\name{MQMLoadanalysis}
\alias{MQMLoadanalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Enables loading of previous MQM analysis from file }
\description{
  Enables loading of previous MQM analysis from file
}
\usage{
MQMLoadanalysis(file = "MQM_output.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ File to load MQM analysis from }
}
\details{

}
\value{
	Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the plot.scanone function of R/QTL.
}
\references{ Danny Arends <Danny.Arends@gmail.com> }
\author{ Danny Arends <Danny.Arends@gmail.com> }
\note{ 	
	If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{  }
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- scanMQM(cross,file="aaa.txt")				# Do a MQM scan of the genome and save it to aaa.txt
res_reloaded <- MQMLoadanalysis("aaa.txt")				#Reads in and plots the file aaa.txt

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
