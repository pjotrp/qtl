\name{CrossFromMolgenis}
\alias{CrossFromMolgenis}

\title{ CrossFromMolgenis - Helperfunction to create a cross-object from a molgenis database }
\description{
  Helperfunction to create a cross-object from a molgenis database
}
\usage{
CrossFromMolgenis(DBmarkerID = 0, DBtraitID = 0, DBpath = NULL)
}

\arguments{
  \item{DBmarkerID}{ Database ID of the matrix with individuals*markers }
  \item{DBtraitID}{ Database ID of the matrix with individuals*traits }
  \item{DBpath}{ Path to the molgenis database }
}
\details{

}
\value{
	Cross object with multiple phenotypes suitable for scanone and scanMQM from R/QTL
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
	\item \code{\link{CrossToMolgenis}} - Function called to store results to a molgenis database
}
}
\examples{
cMol <- CrossFromMolgenis()			#Make a crossobject from the database
result <- scanMQMall(cMol)			#Basic execution of scanMQMall
CrossToMolgenis(result,"QTL-mqm")	#Store the result from MQM into molgenis database
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
