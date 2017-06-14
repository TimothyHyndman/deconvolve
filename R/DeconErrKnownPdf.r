#' Compute the deconvolution KDE
#' 
#' Blahblahblah
#' 
#' PUT DETAILS HERE
#' 
#' @example man/examples/KnownError_eg.R
#' 
#' @export

DeconErrKnownPdf<-function(n,xx,W,h,errortype,sigU,phiU,rescale=0,phiK=phiK2,muK2=6,RK=1024/3003/pi,deltat = .0002,tt = seq(-1,1,deltat))

#Author: Aurore Delaigle
#This function computes the deconvolution kernel density estimator

#----------------------------
#Required arguments:
#----------------------------

#n: sample size
#xx: vector of x-values where to compute the deconvolution kernel density estimator
#W: vector of univariate contaminated data
#h bandwidth

#Error distribution, which can be defined in two ways:

#Option 1: provide error type and standard deviation of the errors
#errortype: 'Lap' for Laplace errors and 'norm' for normal errors. If you use this way of defining the error then you also need to provide the value of sigU below.
#sigU: standard deviation of the errors.

#Option 2: define the characeristic function of the errors:
#phiU:function that gives the characteristic function of the error


#------------------------------------------------------------------------------------------------------------------------
#Optional arguments 
#------------------------------------------------------------------------------------------------------------------------

#rescale: to rescale the estimator so that it integrates to 1 after the negative parts have been truncated to zero. Default: 0= do not rescale. If you want to rescale, set to 1.
#Rescaling requires xx to be a fine grid of equispaced x-values that covers the whole range of x-values where the estimated density is significantly non zero.

#-------------------------------------------------------------------------------------------------------------------------------------
#Changing the kernel: (if you change one of the arguments below you must change them all as they need to correspond to the same kernel):
#-------------------------------------------------------------------------------------------------------------------------------------
#phiK: Fourier transfrom of the kernel. The default is (1-t^2)^3 on the interval [-1,1]
#muK2: second moment of the kernel, i.e. \int x^2 K(x) dx
#RK: integral of the square of the kernel, i.e. \int K^2(x) dx
#tt: vector of discrete t values on which you approximate the integrals in the Fourier domain. 
#	If phiK is compactly supported, the first and last elements of t must be the lower and upper bound of the support of phiK.
#	If phiK is not compactly supported, the first and last elements of t must be larger enough for your discretisation of the intergals to be accurate
#deltat: distance between two points of the t grid 


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#								WARNINGS:
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# The kernel K here must match the phiK used to compute the bandwidth (PI, CV or other)
#
# The DKDE can also be computed using the Fast Fourier Transform, which is a bit more complex. 
# See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution.   Statistics and Computing,  17,  349 - 355
# However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
#
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------




{

#Check optional arguments

if(missing(errortype)&(missing(sigU))&missing(phiU))
	stop("You must define the error distribution")

if(missing(errortype)==F)
{
	if(errortype=='Lap')
		{
		if(missing(sigU))
			stop("You must provide the standard deviation of the errors")

		phiU=function(tt) {return(1/(1+sigU^2*tt^2/2))}
		}

	if(errortype=='norm')
		{
		if(missing(sigU))
			stop("You must provide the standard deviation of the errors")
		phiU=function(tt) {exp(-sigU^2*tt^2/2)}
		}
}


if(missing(h))
	stop("You must provide the bandwidth")


W=as.vector(W)

#make sure t is a vector in the right format
dim(tt)=c(length(tt),1);


OO=outerop(tt/h,t(W),"*")
phiUth=phiU(tt/h)



#Estimate real and imaginary parts of empirical characteristic function of W computed at tt/h, for each component of tt.
rehatphiX=rowSums(cos(OO))/phiUth/n
imhatphiX=rowSums(sin(OO))/phiUth/n

xt=outerop(tt/h,t(xx),"*")
longx=length(xx)

#Compute the DKDE estimator
fXdecUK=cos(xt)*kronecker(matrix(1,1,longx),rehatphiX)+sin(xt)*kronecker(matrix(1,1,longx),imhatphiX)
fXdecUK=apply(fXdecUK*kronecker(matrix(1,1,longx),phiK(tt)),2,sum)/(2*pi)*deltat/h
fXdecUK[which(fXdecUK<0)]=0

if (rescale==1)
	{
	dx=xx[2]-xx[1]
	fXdecUK=fXdecUK/sum(fXdecUK)/dx
	}


return (fXdecUK)}
