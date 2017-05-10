library("deconvolve")
#This code illustrates how to use the functions for computing the deconvolution kernel density estimator and its bandwidths in the case where the errors are heteroscedastic

#-----------------------------------------------------
#Start by generating some data contaminated by noise:
#-----------------------------------------------------


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
tmp=matrix(runif(n,0,1),nrow=1,ncol=n,byrow=TRUE);
X[which(tmp<pmix)]=X2[which(tmp<pmix)];

#Grid where to estimate the true mixture density, and calculation of true density
xx=seq(-5,5,0.1);
dx=xx[2]-xx[1];
truedens=(1-pmix)*dnorm(xx,mu1,sig1)+pmix*dnorm(xx,mu2,sig2);



#Standard devation of the jth error, for j=1,...,n
sigU=0.6;
sigmaj=sigU*sqrt(1.0+(1:n)/n)*sqrt(0.5);

#-------------------------------------
#Example when the error is normal
#------------------------------------
errortype="norm";
U=rep(1,n);
for (i in 1:n)
	U[i]=rnorm(1,0,sigmaj[i]);
 
#Contaminated data
W=as.vector(X+U);


#Estimate the variance of X
varX=mean(W^2)-(mean(W))^2-sum(sigmaj^2)/n;

#PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4het(n,W,varX,errortype,sigmaj);

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknownhet(n,xx,W,hPI,errortype,sigmaj);

#DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknownhet(n,xx,W,hPI,errortype,sigmaj,rescale=1);

#Plot the true density
plot(xx,truedens,'l',col='red',xlab="",ylab="")
lines(xx,y,col='black')
lines(xx,y2,col="green")

#Example of how to provide the vector of phiU_k's instead of the error type and the standard deviations

phiUkvec=c()
for(k in 1:n)
{	
	phiUk<-function(tt,k) {return(exp(-sigmaj[k]^2*tt^2/2));}
	phiUkvec=c(phiUkvec,phiUk)
}

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y3=fdecUknownhet(n,xx,W,hPI,phiUkvec=phiUkvec);
lines(xx,y3,col='magenta',lty=2)



#Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
h=1.06*sqrt(var(W))*n^(-1/5);
xout=outerop(xx,t(W),"-");

fnaive=apply(dnorm(xout,0,h),1,sum)/n;

lines(xx,fnaive,col='cyan')

legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI", "fdec hPI v2", "naive estimator, hNR"),col=c("red","black","green","magenta","cyan"),lty=c(1,1,1,2,1),cex=0.73)