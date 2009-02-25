\name{ResultsFromMolgenis}
\alias{ResultsFromMolgenis}

\title{ ResultsFromMolgenis - Helperfunction to retrieve results from a molgenis database }
\description{
	Helperfunction to retrieve results from a molgenis database
}
\usage{
ResultsFromMolgenis(DBqtlName = NULL, DBpath = NULL)
}

\arguments{
  \item{DBqtlName}{ Name of the matrix holding the QTL results }
  \item{DBpath}{ Path to the molgenis database }
}
\details{

}
\value{
  Returns a MQMmulti or a scanone object.
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
cMol <- CrossFromMolgenis()				#Make a crossobject from the database
result <- scanMQMall(cMol)				#Basic execution of scanMQMall
ResultsToMolgenis(result,"QTL-mqm")			#Store the result from MQM into molgenis database
r_res <- ResultsFromMolgenis("QTL-mqm")		#Store the result from MQM into molgenis database

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }