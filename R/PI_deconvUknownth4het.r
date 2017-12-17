PI_deconvUknownth4het <- function(n, W, varX, phiUkvec, phiK, muK2, RK, deltat,
								  tt){

#Author: Aurore Delaigle
#compute 2-stage plug-in bandwidth for heteroscedastic kerndel deconvolution estimator as in:
# !!! This code is only valid for a kernel of order 2 !!!!

#Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579

#----------------------------
#Required arguments:
#----------------------------

#n: sample size
#W: vector of univariate contaminated data
#varX:estimator of variance of X


#Error distribution, which can be defined in two ways:

#Option 1: provide error type and standard deviation of the errors
#errortype: 'Lap' for Laplace errors and 'norm' for normal errors. If you use this way of defining the error then you also need to provide the value of sigU below.
#sigUj: vector of length n which contains the standard deviations of each of the n errors. 

#Option 2: define the characeristic function of the n errors:
#phiUveck: vector of n functions that give the characteristic functiona of the n errors. Produce this vector by c(func1,func2,...,funcn) where each funcj is a function of tt


#-------------------------------------------------------------------------------------------------------------------------------------
#Changing the kernel: (if you change one of the arguments below you must change them all as they need to correspond to the same kernel):
#-------------------------------------------------------------------------------------------------------------------------------------
# !!! This code is only valid for a kernel of order 2 !!!!
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
# !!! This code is only valid for a kernel of order 2 !!!!
#
# The DKDE can also be computed using the Fast Fourier Transform, which is a bit more complex. 
# See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution.   Statistics and Computing,  17,  349 - 355
# However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
#
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Convert vector of functions to single function -------------------------------
phiUk <- function(tt,k) {
	phiUkvec[[k]](tt)
}


W=as.vector(W)

#grid of h values where to search for a solution: you can change the default grid if no solution is found in this grid.
maxh=(max(W)-min(W))/10

#normal reference bandwidth of the naive KDE estimator (estimator that ignores the errors) using the same kernel as above
hnaive=((8*sqrt(pi)*RK/3/muK2^2)^0.2)*sqrt(stats::var(W))*n^(-1/5)

#grid of h values on which we will look for hPI, If you did not find a minimum on that grid you can redefine it
hgrid=seq(hnaive/3,maxh,(maxh-hnaive/3)/100)
lh = length(hgrid)
dim(hgrid)=c(1,lh)


#Estimator of the standard deviation of X
stdevx = max(sqrt(varX),1/n)

#Quantities that will be needed several times in the computations below
toverh=tt%*%(1/hgrid)

phiKsq=(phiK(tt))^2
phiKsq=as.vector(phiKsq)

#Compute sum over k of phiU_k(tt/h)^2
phiUsq=matrix(0,length(tt),length(hgrid))
for (k in seq_len(n))
	{
	matphiU=phiUk(toverh,k)
	phiUsq=phiUsq+matphiU^2
	}

# --------------------------------------------
# Estimate theta4 by normal reference method     
# --------------------------------------------

th4 = stdevx^(-9)*105/(32*sqrt(pi)) 

# ------------------------------------------------------
# Find bandwidth h3 for computing th3, then compute th3
# ------------------------------------------------------

rr=3
term1= -hgrid^2*muK2*th4
term2=kronecker(matrix(1,1,lh),tt^(2*rr)*phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat
term2=term2/(2*pi*hgrid^(2*rr+1))

ABias2 = (term1 + term2)^2

#Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh3=which.min(ABias2)
if(indh3==1)
	cat("\n minimum of Abias2 for rr=3 is the first element of the grid of bandwidths. Consider enlargin the grid")
if(indh3==length(hgrid))
	cat("\n minimum of Abias2 for rr=3 is the last element of the grid of bandwidths. Consider enlargin the grid")
h3 = hgrid[indh3]

#Estimate empirical characteristic function of W at t/h3

OO <- outer(tt / h3, W)

# Compute phiU(-t/h) -- since phiU is symmetric, this is the same as phiU(t/h)
matphiU=OO
for (k in seq_len(n))
	{matphiU[,k]=phiUk(tt/h3,k)}


phiUsqth=apply(matphiU^2,1,sum)

#Estimate empirical characteristic function of X at t/h3
rehatphiX=apply(cos(OO)*matphiU,1,sum)/phiUsq[,indh3]
imhatphiX=apply(sin(OO)*matphiU,1,sum)/phiUsq[,indh3]

rm(OO)

#Compute th3
normhatphiX2=rehatphiX^2+imhatphiX^2
th3 = sum(tt^(2*rr) * normhatphiX2 * phiKsq)
th3 = th3*deltat/(2*pi*h3^(2*rr+1))

# -----------------------------------------------------
# Find bandwidth h2 for computing th2, then compute th2
# -----------------------------------------------------


rr=2
term1= -hgrid^2*muK2*th3
term2=kronecker(matrix(1,1,lh),tt^(2*rr)*phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat/(2*pi*hgrid^(2*rr+1))
ABias2 = (term1 + term2)^2

#Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh2=which.min(ABias2)
if(indh2==1)
	cat("\n minimum of Abias2 for rr=2 is the first element of the grid of bandwidths. Consider enlargin the grid\n")
if(indh2==length(hgrid))
	cat("\n minimum of Abias2 for rr=2 is the last element of the grid of bandwidths. Consider enlargin the grid\n")
h2 = hgrid[indh2]

#Estimate empirical characteristic function of X at t/h2
OO <- outer(tt / h2, W)
# Compute phiU(-t/h) -- since phiU is symmetric, this is the same as phiU(t/h)
matphiU=OO
for (k in seq_len(n))
	{matphiU[,k]=phiUk(tt/h2,k)}

#Estimate empirical characteristic function of X at t/h2
rehatphiX=apply(cos(OO)*matphiU,1,sum)/phiUsq[,indh2]
imhatphiX=apply(sin(OO)*matphiU,1,sum)/phiUsq[,indh2]

rm(OO)


#Compute th2
normhatphiX2=rehatphiX^2+imhatphiX^2
th2 = sum(tt^(2*rr) * normhatphiX2 * phiKsq)
th2=th2*deltat/(2*pi*h2^(2*rr+1))


# ------------------------------------------------------------------------------------------------------
# Finally, compute the bandwidth that minimises the AMISE of the deconvolution kernel density estimator
# ------------------------------------------------------------------------------------------------------

term1=hgrid^4*muK2^2*th2/4
term2=kronecker(matrix(1,1,lh),phiKsq)/phiUsq
term2=apply(term2,2,sum)*deltat/(2*pi*hgrid)
AMISE=term1+term2

indh=which.min(AMISE)
hPI = hgrid[indh]

return (hPI)}

