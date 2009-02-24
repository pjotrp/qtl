\name{MQMaugment}
\alias{MQMaugment}

\title{ 
	MQMaugment - Missing data augmentation routine for MQM 
}
\description{
  Estimates missing data, unlikely genotypes are kept when more likely than 1:Neglect
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
\references{ Danny Arends \email{Danny.Arends@gmail.com} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{
If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{
	\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
	}
}
\examples{
	library(MQMpackage)
	qtl <- c(3,15,3,7)							# QTL at chromosome 3
	data(map10)								# Mouse genome
	cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
	cross_good <- MQMaugment(cross)						# Augmentation of the data
	result <- scanMQM(cross_good)						# Do a MQM scan of the genome
	plot(result)								# Plot the results of the genome scan
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
