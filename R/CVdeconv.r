CVdeconv<-function(n, W, errortype, sigU, phiU, phiK = phiK2, muK2 = 6, 
				   RK = 1024 / 3003 / pi, deltat = .0002, tt = seq(-1,1,deltat))

# Authors: Aurore Delaigle
# This function computes the cross-validation (CV) bandwidth for kernel 
# deconvolution estimator as in Stefanski, L., Carroll, R.J. (1990). 
# Deconvoluting kernel density estimators. Statistics 2, 169–184. Delaigle, A. 
# and I. Gijbels (2004). Practical bandwidth selection in deconvolution kernel 
# density estimation, Computational Statistics and Data Analysis, 45, 249-267

#----------------------------
# Required arguments:
#----------------------------

# n: sample size
# W: vector of univariate contaminated data

# Error distribution, which can be defined in two ways:

# Option 1: provide error type and standard deviation of the errors
# errortype: 'Lap' for Laplace errors and 'norm' for normal errors. If you use 
# this way of defining the error then you also need to provide the value of 
# sigU below.
# sigU: standard deviation of the errors.

# Option 2: define the characeristic function of the errors:
# phiU:function that gives the characteristic function of the error


#------------------------------------------------------------------------------------------------------------------------
# Optional arguments (if you change one of them you must change all of them as they need to correspond to the same kernel):
#------------------------------------------------------------------------------------------------------------------------


#phiK: Fourier transfrom of the kernel. The default is (1-t^2)^3 on the interval [-1,1]
#muK2: second moment of the kernel, i.e. \int x^2 K(x) dx
#RK: integral of the square of the kernel, i.e. \int K^2(x) dx
#tt: vector of discrete t values on which you approximate the integrals in the Fourier domain. 
#	If phiK is compactly supported, the first and last elements of t must be the lower and upper bound of the support of phiK.
#	If phiK is not compactly supported, the first and last elements of t must be larger enough for your discretisation of the intergals to be accurate
#deltat: distance between two points of the t grid 


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#								WARNINGS:
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#The kernel you use must be the same as the kernel defined in the function fdecUknown, so if you change the optional arguments here you must change them in fdecUknown.
#
#If you change the kernel you have to chage muK2, RK and the range of t-values (these must correspond to the support of phiK), and thus also delatat
#
#In case of multiple bandwidth solutions, by default this code takes the largest solution: you can change this to your preferred way of breaking ties.
#Often if you plot CV you will see that the first few solutions seem unreasonable (CV fluctuates widely). You can take the first minimum that looks reasonable.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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



# --------------------------------------------------------
# Preliminary calculations and initialisation of functions
# --------------------------------------------------------


W=as.vector(W)

#make sure t is a vector in the right format
dim(tt)=c(length(tt),1);


#Define hgrid, the grid of h values where to search for a solution: you can change the default grid if no solution is found on this grid.
maxh=(max(W)-min(W))/10;

#normal reference bandwidth of the naive KDE estimator (estimator that ignores the errors) using the same kernel as above
hnaive=((8*sqrt(pi)*RK/3/muK2^2)^0.2)*sqrt( stats::var(W) )*n^(-1/5);

#grid of h values on which we will look for hCV, If you did not find a minimum on that grid you can redefine it
hgrid=seq(hnaive/3,maxh,(maxh-hnaive/3)/100);
lh = length(hgrid);
dim(hgrid)=c(1,lh);


#Quantities that will be needed several times in the computations below
toverh=tt%*%(1/hgrid);
phiU2=phiU(toverh)^2;
phiKt=phiK(tt);



#----------------------
#Compute CV criterion
#----------------------

longh=length(hgrid);
CVcrit=rep(0,longh);
OO=outerop(tt,t(W),"*");


#Compute CV criterion for all values of h on the grid of h values

for (j in 1:longh)
	{
	h=hgrid[j];

	#Estimate the square of the norm of the empirical characteristic function of W
	rehatphiW=apply(cos(OO/h),1,sum)/n;
	imhatphiW=apply(sin(OO/h),1,sum)/n;
	normhatphiW2=rehatphiW^2+imhatphiW^2;
	
	#Compute CV
	CVcrit[j]=sum(phiKt/phiU2[,j]*(normhatphiW2*((n-1)*phiKt-2*n)+2));
	CVcrit[j]=CVcrit[j]/h;
	}


#----------------------
#Find CV bandwidth
#----------------------

#Find the indices corresponding to all local minima of the CV curve
indh=which((CVcrit[2:(longh-1)]<CVcrit[1:(longh-2)])&(CVcrit[2:(longh-1)]<CVcrit[3:(longh)]));


#if did not find a local minimum, take the global minimun on the boundaries of the h grid
if (length(indh)<1)
	{hCV=which.min(CVcrit)} else {
	#In case of multiple solutions, take the largest bandwidth: you can change this to your preferred way of breaking ties
	hCV = max(hgrid[indh]);}


return (hCV)}

