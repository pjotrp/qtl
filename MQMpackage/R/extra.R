extra <- function(a,b,x){
	.C("R_gammaln",as.double(3.0))
	lgamma(a)

	.C("R_betacf",as.double(3.0),as.double(3.0),as.double(0.9))

	.C("r_betai",a= as.double(3.0),b= as.double(3.0),x= as.double(0.9))
	pbeta(x,a,b)
}