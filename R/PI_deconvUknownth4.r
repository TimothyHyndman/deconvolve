#' @export

PI_deconvUknownth4<-function(n, W, errortype, sigU, phiU, phiK = phiK2, 
							 muK2 = 6, RK = 1024 / 3003 / pi, deltat = .0002, 
							 tt = seq(-1, 1, deltat)){

# Authors: Aurore Delaigle
# compute 2-stage plug-in bandwidth for kerndel deconvolution estimator as in:
# Delaigle, A. and I. Gijbels (2002). Estimation of integrated squared density 
# derivatives from a contaminated sample, Journal of the Royal Statistical 
# Society, B, 64, 869-886.
# Delaigle, A. and I. Gijbels (2004). Practical bandwidth selection in 
# deconvolution kernel density estimation, Computational Statistics and Data 
# Analysis, 45, 249 - 267

# !!! This code is only valid for a kernel of order 2 !!!!


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


#-------------------------------------------------------------------------------
# Optional arguments (if you change one of them you must change all of them as 
# they need to correspond to the same kernel):
#-------------------------------------------------------------------------------

# !!! This code is only valid for a kernel of order 2 !!!!
# phiK: Fourier transfrom of the kernel. The default is (1-t^2)^3 on the 
# interval [-1,1]
# muK2: second moment of the kernel, i.e. \int x^2 K(x) dx
# RK: integral of the square of the kernel, i.e. \int K^2(x) dx
# tt: vector of discrete t values on which you approximate the integrals in the 
# Fourier domain. 
#	If phiK is compactly supported, the first and last elements of t must be the 
# 	lower and upper bound of the support of phiK.
#	If phiK is not compactly supported, the first and last elements of t must be 
#	larger enough for your discretisation of the intergals to be accurate
# deltat: distance between two points of the t grid 


#-------------------------------------------------------------------------------
#								WARNINGS:
#-------------------------------------------------------------------------------
# !!! This code is only valid for a kernel of order 2 !!!!
# The kernel you use must be the same as the kernel defined in the function 
# fdecUknown, so if you change the optional arguments here you must change them 
# in fdecUknown.
#
# If you change the kernel you have to chage muK2, RK and the range of t-values 
# (these must correspond to the support of phiK), and thus also delatat
#
# In case of multiple bandwidth solutions, by default this code takes the 
# largest solution: you can change this to your preferred way of breaking ties.
# Often if you plot CV you will see that the first few solutions seem 
# unreasonable (CV fluctuates widely). You can take the first minimum that looks 
# reasonable.
#-------------------------------------------------------------------------------


# Check optional arguments

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


# second moment \int x^2 K(x) dx  of the kernel K
W=as.vector(W)


# grid of h values where to search for a solution: you can change the default 
# grid if no solution is found in this grid.
maxh=(max(W)-min(W))/10

# normal reference bandwidth of the naive KDE estimator (estimator that ignores 
# the errors) using the same kernel as above
hnaive=((8*sqrt(pi)*RK/3/muK2^2)^0.2)*sqrt(stats::var(W))*n^(-1/5)

# grid of h values on which we will look for hPI, If you did not find a minimum 
# on that grid you can redefine it
hgrid=seq(hnaive/3,maxh,(maxh-hnaive/3)/100)
lh = length(hgrid)
dim(hgrid)=c(1,lh)


# Estimator of the standard deviation of X
stdevx = max(sqrt(stats::var(W) - sigU^2),1/n) ##### I changed varU into sigU^2



# Quantities that will be needed several times in the computations below
toverh=tt%*%(1/hgrid)

phiKsq=(phiK(tt))^2
phiKsq=as.vector(phiKsq)
phiUsq=phiU(toverh)^2




#---------------------------------------------
# Estimate theta4 by normal reference method     
#---------------------------------------------

th4 = stdevx^(-9)*105/(32*sqrt(pi)) 

#-------------------------------------------------------
# Find bandwidth h3 for computing th3, then compute th3
#-------------------------------------------------------

rr=3
term1= -hgrid^2*muK2*th4
term2=kronecker(matrix(1,1,lh),tt^(2*rr)*phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat
term2=term2/(2*pi*n*hgrid^(2*rr+1))

ABias2 = (term1 + term2)^2

# Print the index of the minimiser of Abias2 to see if we are inside the grid 
# (if not, enlarge the grid of bandwidths)
indh3=which.min(ABias2)
if(indh3==1)
	cat("\n minimum of Abias2 for rr=3 is the first element of the grid of 
		bandwidths. Consider enlarging the grid")
if(indh3==length(hgrid))
	cat("\n minimum of Abias2 for rr=3 is the last element of the grid of 
		bandwidths. Consider enlarging the grid")
h3 = hgrid[indh3]

# Estimate empirical characteristic function of W at t/h3
OO=outerop(tt/h3,t(W),"*")
rehatphiW=apply(cos(OO),1,sum)/n
imhatphiW=apply(sin(OO),1,sum)/n
rm(OO)

# Compute th3
normhatphiW2=rehatphiW^2+imhatphiW^2
th3 = sum(tt^(2*rr) * normhatphiW2 * phiKsq / phiUsq[,indh3])
th3 = th3*deltat/(2*pi*h3^(2*rr+1))

# -----------------------------------------------------
# Find bandwidth h2 for computing th2, then compute th2
# -----------------------------------------------------


rr=2
term1= -hgrid^2*muK2*th3
term2=kronecker(matrix(1,1,lh),tt^(2*rr)*phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat/(2*pi*n*hgrid^(2*rr+1))

ABias2 = (term1 + term2)^2

# Print the index of the minimiser of Abias2 to see if we are inside the grid 
# (if not, enlarge the grid of bandwidths)
indh2=which.min(ABias2)
if(indh2==1)
	cat("\n minimum of Abias2 for rr=2 is the first element of the grid of 
		bandwidths. Consider enlarging the grid\n")
if(indh2==length(hgrid))
	cat("\n minimum of Abias2 for rr=2 is the last element of the grid of 
		bandwidths. Consider enlarging the grid\n")
h2 = hgrid[indh2]

# Estimate empirical characteristic function of W at t/h2
OO=outerop(tt/h2,t(W),"*")
rehatphiW=apply(cos(OO),1,sum)/n
imhatphiW=apply(sin(OO),1,sum)/n
rm(OO)


# Compute th2
normhatphiW2=rehatphiW^2+imhatphiW^2
th2 = sum(tt^(2*rr) * normhatphiW2 * phiKsq / phiUsq[,indh2])
th2=th2*deltat/(2*pi*h2^(2*rr+1))


#-------------------------------------------------------------------------------
# Finally, compute the bandwidth that minimises the AMISE of the deconvolution 
# kernel density estimator
#-------------------------------------------------------------------------------

term1=hgrid^4*muK2^2*th2/4


term2=kronecker(matrix(1,1,lh),phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat/(2*pi*n*hgrid)
AMISE=term1+term2

indh=which.min(AMISE)
hPI = hgrid[indh]

hPI
}