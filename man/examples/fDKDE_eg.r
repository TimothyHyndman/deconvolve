#Put here the path where you have put the R codes, for example:

library(fDKDE)

source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/CVdeconv.R")
source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/PI_deconvUknownth4.r")
source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/fdecUknown.r")
source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/phiK2.r")
source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/rlap.r")
source("C:/Users/delaigle/Documents/articles post these/CodesDeconvolution/decon_Rpackages_Tianying/code_decon/fDKDE/outerop.R")


#Author: Aurore Delaigle
#This code illustrates how to use the functions for computing the deconvolution kernel density estimator and its bandwidths

#-----------------------------------------------------
#Start by generating some data contaminated by noise:
#-----------------------------------------------------



#Noise to signal ratio=varU/varX
NSR=0.2

#Sample size
n=200

#Generate data from a normal mixture
mu1=-3;
mu2=2;
sig1=1;
sig2=1;

X=rnorm(n,mu1,sig1);
X2=rnorm(n,mu2,sig2);

pmix=0.75;
tmp=matrix(runif(n,0,1),nrow=1,ncol=n,byrow=T);
X[which(tmp<pmix)]=X2[which(tmp<pmix)];

#Grid where to estimate the true mixture density, and calculation of true density
xx=seq(-5,5,0.1);
dx=xx[2]-xx[1];
truedens=(1-pmix)*dnorm(xx,mu1,sig1)+pmix*dnorm(xx,mu2,sig2);



#-------------------------------------
#Example when the error is normal
#-------------------------------------

errortype="norm";
sigU=sqrt(NSR*var(X));
U=rnorm(n,0,sigU);
W=as.vector(X+U);

#Plot the true density
plot(xx,truedens,'l',col='red',xlab="",ylab="")

#PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4(n,W,errortype,sigU);

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknown(n,xx,W,hPI,errortype,sigU);

#DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknown(n,xx,W,hPI,errortype,sigU,rescale=1);

lines(xx,y2,col="green",xlab="",ylab="")
lines(xx,y,col='black')



#CV bandwidth of Stefanski and Carroll
hCV=CVdeconv(n,W,errortype,sigU)

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y3=fdecUknown(n,xx,W,hCV,errortype,sigU);

lines(xx,y3,col='magenta')



#Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
h=1.06*sqrt(var(W))*n^(-1/5);
xout=outerop(xx,t(W),"-");

fnaive=apply(dnorm(xout,0,h),1,sum)/n;

lines(xx,fnaive,col='cyan')


legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI", "fdec rescaled, hCV", "naive estimator, hNR"),col=c("red","black","green","magenta","cyan"),lty=c(1,1,1,1),cex=0.73)



#-------------------------------------
#Example when the error is Laplace
#-------------------------------------
windows()
errortype="Lap"
sigLap=sqrt(NSR*var(X)/2)
sigU=sqrt(2)*sigLap;
U=rlap(sigLap,1,n);

#Contaminated data
W=as.vector(X+U);


#Plot the true density
plot(xx,truedens,'l',col='red',xlab="",ylab="")

#PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4(n,W,errortype,sigU);


#DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknown(n,xx,W,hPI,errortype,sigU);

#DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknown(n,xx,W,hPI,errortype,sigU,rescale=1);

lines(xx,y2,col="green",xlab="",ylab="")
lines(xx,y,col='black')



#CV bandwidth of Stefanski and Carroll
hCV=CVdeconv(n,W,errortype,sigU)

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y3=fdecUknown(n,xx,W,hCV,errortype,sigU);

lines(xx,y3,col='magenta')



#Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
h=1.06*sqrt(var(W))*n^(-1/5);
xout=outerop(xx,t(W),"-");

fnaive=apply(dnorm(xout,0,h),1,sum)/n;

lines(xx,fnaive,col='cyan')


legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI", "fdec rescaled, hCV", "naive estimator, hNR"),col=c("red","black","green","magenta","cyan"),lty=c(1,1,1,1),cex=0.73)



