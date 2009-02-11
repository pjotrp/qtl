\name{MQMlogPheno}
\alias{MQMlogPheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MQMlogPheno(cross = NULL, Phenot = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ ~~Describe \code{cross} here~~ }
  \item{Phenot}{ ~~Describe \code{Phenot} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(cross= NULL,Phenot=NULL){
    library(qtl)
if(is.null(cross)){
stop("Error: No cross file. Please supply a valid cross object.")
return 
  }
if(class(cross)[1] == "f2" || class(cross)[1] == "bc" || class(cross)[1] == "ril"){
if(!is.null(Phenot)){
pheno <- NULL
logpheno <- NULL
pheno <- cross$pheno[[Phenot]]
logpheno <- log(pheno)
cross$pheno[[Phenot]] <- logpheno
cat("Phenotype:",Phenot,".\n")
cat("Before LOG MEAN:",mean(pheno,na.rm = TRUE),"VAR:",var(pheno,na.rm = TRUE),".\n")
cat("After LOG MEAN:",mean(logpheno,na.rm = TRUE),"VAR:",var(logpheno,na.rm = TRUE),".\n")
  }else{
n.pheno <- nphe(cross)
for(i in 1:n.pheno) {
pheno <- cross$pheno[[i]]
logpheno <- log(pheno)
cross$pheno[[i]] <- logpheno
cat("Phenotype:",i,".\n")
cat("Before LOG MEAN:",mean(pheno,na.rm = TRUE),"VAR:",var(pheno,na.rm = TRUE),".\n")
cat("After LOG MEAN:",mean(logpheno,na.rm = TRUE),"VAR:",var(logpheno,na.rm = TRUE),".\n")
  }
  }
cross
  }else{
stop("Error: Currently only F2 / BC / RIL cross files can be analyzed by MQM.")
  }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
