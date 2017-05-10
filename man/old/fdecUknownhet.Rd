\name{fdecUknownhet}
\alias{fdecUknownhet}
\title{
	Compute the deconvolution kernel density estimator
}
\description{
	Compute the deconvolution kernel density estimator when the errors are heteroscedastic, as in Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579
}
\usage{
	fdecUknownhet(n,xx,W,h,errortype,sigUj,phiUkvec,rescale=0,phiK=phiK2,
	muK2=6,RK=1024/3003/pi,deltat = .0002,tt = seq(-1,1,deltat))
}
\arguments{
	\item{n}{
	Sample size
	}
	\item{xx}{
	vector of x-values where to compute the deconvolution kernel density estimator
	}
	\item{W}{
	vector of univariate contaminated data
	}
	\item{h}{
	bandwidth
	}
	\item{errortype}{
	"Lap" for the case where the error densities are all Laplace densities and "norm" for the case where the error densities are all normal. If you use this way of defining the error then you also need to provide the value of sigUj below.
	}
	\item{sigUj}{
	vector of length n which contains the parameters of each of the n Laplace or normal errors.
	}
	\item{phiUveck}{
	vector of n functions that give the characteristic functiona of the n errors. Produce this vector by c(func1,func2,...,funcn) where each funcj is a function of tt
	}
	\item{rescale}{
	to rescale the estimator so that it integrates to 1 after the negative parts have been truncated to zero. Default is 0 (do not rescale). If you want to rescale, set to 1 and see more details below.
	}
	\item{phiK}{
	Fourier transfrom of the kernel. The default is \eqn{(1-t^2)^3} on the interval \eqn{[-1,1]}
	}
	\item{muK2}{
	second moment of the kernel, i.e. \eqn{x^2 K(x) dx}
	}
	\item{RK}{
	integral of the square of the kernel, i.e. \eqn{ K^2(x) dx}
	}
	\item{deltat}{
	distance between two points of the t grid.
	}
	\item{tt}{
	vector of discrete t values on which you approximate the integrals in the Fourier domain.
	}
}
\section{Warnings}{
	\enumerate{
	\item Rescaling requires xx to be a fine grid of equispaced x-values that covers the whole range of x-values where the estimated density is significantly non zero.
	\item Changing the kernel: if you change one of the arguments among phiK, muK2, RK, deltat and tt, you must change them all as they need to correspond to the same kernel.
	\item	If phiK is compactly supported, the first and last elements of t must be the lower and upper bound of the support of phiK.
	\item	If phiK is not compactly supported, the first and last elements of t must be larger enough for your discretisation of the intergals to be accurate
	\item The kernel K here must match the phiK used to compute the bandwidth (PI, CV or other)
	\item The DKDE can also be computed using the Fast Fourier Transform, which is a bit more complex. See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution. Statistics and Computing, 17, 349-355
	\item However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
	}
}
\value{The outcome is the deconvolution kernel density estimator when the errors are heteroscedastic.
}
\references{
	Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579}
\author{
	Aurore Delaigle}
\examples{
#See examples for this package
}
\keyword{deconvolution kernel density estimator}
\keyword{heteroscedastic}
