\name{MQMpackage-package}
\alias{MQMpackage-package}
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
MQMaugment              MQMaugment - Missing data augmentation routine
                        for MQM
MQMCofactors            MQMCofactors - Helper function to create
                        cofactors to be used with scanMQM
MQMCofactorsEach        MQMCofactorsEach - Helper function to create
                        cofactors to be used with scanMQM
multitrait              Crossobject from R/QTL with multiple traits
scanMQM                 scanMQM - Multiple qtl mapping for F2 crosses
scanMQMall              scanMQMall - Multitrait analysis of a cross
                        object
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Danny Arends <Danny.Arends@gmail.com> ; R.C. Janssen <> ; K. Brohamn <>

Maintainer: Danny Arends <Danny.Arends@gmail.com>
}
\references{

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{

}
\examples{
	library(MQMpackage)
	data(multitrait)
	result <- scanMQM(multitrait)
	plot(result)
}
