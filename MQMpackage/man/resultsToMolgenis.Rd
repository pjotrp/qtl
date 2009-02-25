\name{ResultsToMolgenis}
\alias{ResultsToMolgenis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ResultsToMolgenis - Helperfunction to store results from a MQM-qtl run into a molgenis database }
\description{
  Helperfunction to store results from a MQM-qtl run into a a molgenis database
}
\usage{
ResultsToMolgenis(intervalQTLmap = NULL, name = NULL, DBpath = NULL, Fupdate = 0)
}
\arguments{
  \item{intervalQTLmap}{ Results from either scanone, scanMQM or scanMQMall }
  \item{name}{ Name of the table to store the results in }
  \item{DBpath}{ Path to the molgenis database }
  \item{Fupdate}{ Force update of values in the database ? (0= don't update / 1= update) }
}
\details{
  
}
\value{
	No return value
}
\references{ http://www.molgenis.org \& http://www.Xgap.org }
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; M.A. Swertz \email{m.a.swertz@rug.nl} }

\note{
The library RCurl from bioconductor should be installed.
\itemize{
	\item http://www.molgenis.org
	\item http://www.Xgap.org
}
If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do single trait analysis
	\item \code{\link{CrossFromMolgenis}} - Function called to store results to a molgenis database
}
}
\examples{
cMol <- CrossFromMolgenis()			#Make a crossobject from the database
result <- scanMQMall(cMol)			#Basic execution of scanMQMall
ResultsToMolgenis(result,"QTL-mqm")	#Store the result from MQM into molgenis database
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
