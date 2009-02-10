\name{plot.MQMall}
\alias{plot.MQMall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plot.MQMall - Plotting routine to display the results from a scanMQMall scan }
\description{
  Plotting routine to display the results from a scanMQMall scan, in either a contourmap or a multiple QTL plot from R/QTL plot.scanone
}
\usage{
plot.MQMall(cross = NULL, result = NULL, type = "C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{result}{ Result object from scanMQMall }
  \item{type}{ either "C" for contourplot or "P" for stacked R/QTL plot }
}
\details{
  
}
\value{
	No return, plotting routine
}
\references{  }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com)
}
\seealso{ 
	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
\examples{
	library(MQMpackage)
	data(multitrait)
	result <- scanMQMall(multitrait)
	plot.MQMall(cross,result,"P")
	plot.MQMall(cross,result,"C")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
