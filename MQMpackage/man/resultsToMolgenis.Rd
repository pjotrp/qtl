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
 The library RCurl from bioconductor should be installed to use a molgenis interface.
}
\value{
	No return value
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; M.A. Swertz \email{m.a.swertz@rug.nl} ; R.C. Janssen \email{r.c.jansen@rug.nl} }

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
\references{
	http://www.molgenis.org \& http://www.Xgap.org	
	\itemize{
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Swertz MA, De Brock EO, Van Hijum SA, De Jong A, Buist G, Baerends RJ, Kok J, Kuipers OP, Jansen RC. (2004), Molecular Genetics Information System (MOLGENIS): alternatives in developing local experimental genomics databases. \emph{Bioinformatics},\bold{13}, 2075--2083.
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
