\name{MQMpackage-package}
\alias{MQMpackage}
\docType{package}
\title{
Tool for analyzing QTL experiments
}
\description{
Analysis of experimental crosses to identify QTLs contributing to variation in quantitative traits.
}
\details{
\tabular{ll}{
Package: \tab MQMpackage\cr
Version: \tab 1.0\cr
Date: \tab 30 January 2009\cr
Depends: \tab R (>= 1.8.0), qtl\cr
License: \tab GPL version 2 or later\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Packaged: \tab Mon Feb  9 13:36:40 2009; F101715\cr
Built: \tab R 2.8.1; i386-pc-mingw32; 2009-02-09 13:36:42; windows\cr
}

Index:
\preformatted{
scanMQM                 scanMQM - Multiple QTL mapping for F2 crosses
scanMQMall              scanMQMall - Multitrait analysis of a cross object
bootstrapMQM		bootstrapMQM - Bootstrapping approach to obtain estimates of QTL significance
MQMaugment              MQMaugment - Missing data augmentation routine for MQM
MQMCofactors            MQMCofactors - Helper function to create cofactors to be used with scanMQM
MQMCofactorsEach        MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM
CrossFromMolgenis       CrossFromMolgenis - Helperfunction to create a cross-object from a molgenis database
ResultsFromMolgenis     ResultsFromMolgenis - Helperfunction to retrieve results from a molgenis database
ResultsToMolgenis       ResultsToMolgenis - Helperfunction to store results into a molgenis database
multitrait              Crossobject from R/QTL with multiple traits
}
Normal layout of MQM-QTL analysis
\itemize{
\item Load in experiment using formats supported by R/QTL
\item Set trait under consideration as first trait in R/QTL format
\item Fill in missing data using either the fill.geno or MQMaugmentdata routine
\item Call scanMQM to do a genomescan
\item Call scanMQM setting the n.run to get an estimate of ‘significant’ lod scores
\item Set cofactors at ‘top scoring’ and not so ‘top scoring’ markers
\item Rerun analysis with cofactors
\item Bootstrap analysis with cofactors to get estimates of QTL significance
}
}
\author{
Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} ; M.A. Swertz \email{m.a.swertz@rug.nl}

Maintainer: Danny Arends \email{Danny.Arends@gmail.com}
}
\references{

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }

\seealso{
	\itemize{
	\item	\code{\link{MQMaugment}} - Dataaugmentation routine to estimate missing data
	\item	\code{\link{scanMQM}} - Function called to do single trait analysis
	\item	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
	}
}
\examples{
	library(MQMpackage)
	data(multitrait)
	mt_est <- MQMaugment(multitrait)
	result <- scanMQM(mt_est)
	plot(result)
}
