\name{MQMfind.marker}
\alias{MQMfind.marker}
\title{ Find the genetic markers significant after bootstrapping/permutation of the data }
\description{
 Find the genetic markers significant after bootstrapping/permutation of the data
}
\usage{
MQMfind.marker(cross = NULL, scanMQMres = NULL, perm = NULL, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and scanMQM from R/QTL. See \code{\link[qtl]{read.cross}} for details.  }
  \item{scanMQMres}{ Results from either scanone or scanMQM }
  \item{perm}{ a scanoneperm object }
  \item{alpha}{ Threshold value, everything with significance < alpha is reported }
}
\details{

}
\value{
	return a matrix with at each row a significant marker (determined from the scanoneperm object) and with columns: markername, chr and pos (Cm)
}
\references{ 

}
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link[qtl]{find.marker}} - R/QTL find.marker routine
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
	\item \code{\link{MQMpermObject}} - Function called to convert MQMbootstrap into an scanoneperm object
}
}
\examples{
library(MQMpackage)
f2qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")		# Simulate a F2 Cross
f2result <- bootstrapMQM(f2cross)				# Do a MQM scan of the genome
f2perm <- MQMpermObject(f2result)				#Create a permutation object
summary(f2perm)				#What LOD score is considered significant ?
numresults <- dim(f2result[[1]])[1]
lines(x=1:numresults,y=rep(summary(f2perm)[1,1],numresults),col="green",lwd=2,lty=2)	#Green line for A=0.05
lines(x=1:numresults,y=rep(summary(f2perm)[2,1],numresults),col="blue",lwd=2,lty=2)		#Blue line for A=0.10
marker <- MQMfind.marker(f2cross,f2result[[1]],f2perm)		#Find markers with a significant QTL effect (First run is original phenotype data)
marker								#Print it to the screen
lodint(f2result[[1]],chr=marker[,'chr']) 			#Extra: Confidence intervals
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
