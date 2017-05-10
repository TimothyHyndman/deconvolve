#' Compute the deconvolution kernel density estimator.
#' 
#' Compute the deconvolution kernel density estimator when the errors are 
#' heteroscedastic, as in Delaigle, A. and Meister, A. (2008). Density 
#' estimation with heteroscedastic error. Bernoulli, 14, 562-579
#' 
#' PUT DETAILS HERE
#' 
#' @param n Sample size
#' @param xx vector of x-values where to compute the deconvolution kernel 
#' density estimator
#' @param W vector of univariate contaminated data
#' @param h bandwidth
#' @param errortype "Lap" for the case where the error densities are all Laplace
#' densities and "norm" for the case where the error densities are all normal. 
#' If you use this way of defining the error then you also need to provide the 
#' value of sigUj below.
#' @param sigUj	vector of length n which contains the parameters of each of the 
#' n Laplace or normal errors.
#' @param phiUveck vector of n functions that give the characteristic functiona 
#' of the n errors. Produce this vector by c(func1,func2,...,funcn) where each 
#' funcj is a function of tt
#' @param rescale to rescale the estimator so that it integrates to 1 after the 
#' negative parts have been truncated to zero. Default is 0 (do not rescale). 
#' If you want to rescale, set to 1 and see more details below.
#' @param phiK Fourier transfrom of the kernel. The default is \eqn{(1-t^2)^3} 
#' on the interval \eqn{[-1,1]}
#' @param muK2 second moment of the kernel, i.e. \eqn{x^2 K(x) dx}
#' @param RK integral of the square of the kernel, i.e. \eqn{ K^2(x) dx}
#' @param deltat distance between two points of the t grid.
#' @param tt vector of discrete t values on which you approximate the integrals 
#' in the Fourier domain.
#' 
#' @section Warnings:
#' \enumerate{
#'	\item Rescaling requires xx to be a fine grid of equispaced x-values that 
#'	covers the whole range of x-values where the estimated density is 
#'	significantly non zero.
#'	\item Changing the kernel: if you change one of the arguments among phiK, 
#'	muK2, RK, deltat and tt, you must change them all as they need to correspond 
#'	to the same kernel.
#'	\item	If phiK is compactly supported, the first and last elements of t 
#'	must be the lower and upper bound of the support of phiK.
#'	\item	If phiK is not compactly supported, the first and last elements of t 
#'	must be larger enough for your discretisation of the intergals to be 
#'	accurate
#'	\item The kernel K here must match the phiK used to compute the bandwidth 
#'	(PI, CV or other)
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle, A. and Gijbels, I. (2007). Frequent 
#' 	problems in calculating integrals and optimizing objective functions: a case 
#' 	study in density deconvolution. Statistics and Computing, 17, 349-355
#'	\item However if the grid of t-values is fine enough, the estimator can 
#' 	simply be computed like here without having problems with oscillations.
#'	}
#' 
#' @return The outcome is the deconvolution kernel density estimator when the 
#' errors are heteroscedastic.
#' 
#' @section References:
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579
#' 
#' @section Author:
#' Aurore Delaigle
#' 
#' @examples #See examples for this package
#' 
#' @keywords heteroscedastic deconvolution kernel density estimator

fdecUknownhet<-function(n,xx,W,h,errortype,sigUj,phiUkvec,rescale=0,phiK=phiK2,muK2=6,RK=1024/3003/pi,deltat = .0002,tt = seq(-1,1,deltat))


#Author: Aurore Delaigle
#Compute the deconvolution kernel density estimator when the errors are heteroscedastic.
#as in Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579

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
#sigUj: vector of length n which contains the standard deviations of each of the n errors. 

#Option 2: define the characeristic function of the n errors:
#phiUveck: vector of n functions that give the characteristic functiona of the n errors. Produce this vector by c(func1,func2,...,funcn) where each funcj is a function of tt


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

if(missing(errortype)&missing(phiUkvec))
	stop("You must define the error distributions")

if(missing(errortype)==F)
{
	if(errortype=='Lap')
		{
		if(missing(sigUj))
			stop("You must provide the standard deviations of the errors")

		phiUk=function(tt,k) {return(1/(1+sigUj[k]^2*tt^2/2))}
		}

	if(errortype=='norm')
		{
		if(missing(sigUj))
			stop("You must provide the standard deviations of the errors")
		phiUk=function(tt,k) {exp(-sigUj[k]^2*tt^2/2)}
		}
}


#characteristic function of kth error 
if(missing(errortype)&(missing(phiUkvec)==F))
		phiUk=function(tt,k)  {phiUkvec[[k]](tt,k)}


	

	
	
if(missing(h))
	stop("You must provide the bandwidth")


W=as.vector(W)

#make sure t is a vector in the right format
dim(tt)=c(length(tt),1);






#Default values of phiU(t)=characteristic function of the errors
#If you want to consider another error type, simply replace phiU by the characteristic function of your error type


OO=outerop(tt/h,t(W),"*")

# Compute phiU_k(-t/h) for each k -- since phiU_k is symmetric, this is the same as phiU_k(t/h)
matphiU=OO;
for (k in 1:n)
	matphiU[,k]=phiUk(tt/h,k)

#sum by rows of the matrix. This produces a vector of size equal to that of tt
phiUsqth=apply(matphiU^2,1,sum)



#Estimate real and imaginary parts of empirical characteristic function of W computed at tt/h, for each component of tt.
#Results=vectors of size length(tt)

rehatphiX=apply(cos(OO)*matphiU,1,sum)/phiUsqth
imhatphiX=apply(sin(OO)*matphiU,1,sum)/phiUsqth


#Matrix of size length(tt) x length(xx)
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

return (fXdecUK)

}
