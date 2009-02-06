\name{MQMaugment}
\alias{MQMaugment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MQMaugment - Missing data augmentation routine for MQM }
\description{
  
}
\usage{
MQMaugment(cross = NULL, maxaug = 10000, maxiaug = 1000, neglect = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{Phenot}{ Phenotype to pull from the cross object }  
  \item{maxaug}{ Maximum number of (total) individuals after augmentation }
  \item{maxiaug}{ Maximum augmentations per individual }
  \item{neglect}{ Genotypes considered have a higher probability of occuring then one in neglect }
}
\details{
Special care has 2 be taken when setting unrealistic values.
When the algorithm ends with an "DATASET TOO LARGE after augmentation", the crossfile could be corrupted.
}
\value{
Cross object with 1 phenotype (indicated by setting the Phenot parameter) suitable for scanMQM from R/QTL
}
\references{ ~put references to the literature/web site here ~ }
\author{ Danny Arends }
\note{
If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	library(MQMpackage)
	qtl <- c(3,15,3,7)							# QTL at chromosome 3
	data(map10)								# Mouse genome
	cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
	cross_good <- MQMaugment(cross)						# Augmentation of the data
	result <- scanMQM(cross_good)						# Do a MQM scan of the genome
	plot(result)								# Plot the results of the genome scan
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
