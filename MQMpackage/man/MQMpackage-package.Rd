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
MQMaugment              MQMaugment - Missing data augmentation routine for MQM
MQMCofactors            MQMCofactors - Helper function to create cofactors to be used with scanMQM
MQMCofactorsEach        MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM
multitrait              Crossobject from R/QTL with multiple traits
scanMQM                 scanMQM - Multiple QTL mapping for F2 crosses
scanMQMall              scanMQMall - Multitrait analysis of a cross object
}
Normal layout of MQM-QTL analysis
-Load in experiment using formats supported by R/QTL
-Set trait under consideration as first trait in R/QTL format
-Fill in missing data using either the fill.geno or MQMaugmentdata routine
-Call scanMQM to do a genomescan
-Call scanMQM setting the n.run to get an estimate of ‘significant’ lod scores
-Set cofactors at ‘top scoring’ and not so ‘top scoring’ markers
-Rerun analysis with cofactors

}
\author{
Danny Arends <Danny.Arends@gmail.com> ; R.C. Janssen <r.c.jansen@rug.nl> ; Karl W. Brohamn <kbroman@biostat.wisc.edu>

Maintainer: Danny Arends <Danny.Arends@gmail.com>
}
\references{

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
	\code{\link{MQMaugment}} - Dataaugmentation routine to estimate missing data
	\code{\link{scanMQM}} - Function called to do single trait analysis
	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
\examples{
	library(MQMpackage)
	data(multitrait)
	mt_est <- MQMaugment(multitrait)
	result <- scanMQM(mt_est)
	plot(result)
}
