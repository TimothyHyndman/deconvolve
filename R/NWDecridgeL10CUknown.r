##This function is seperated from hSIMEX.r
NWDecridgeL1OCUknown<-function(n,W,Y,errortype,sigU,h,rhogrid,midbin,indbin,nbin)
{
  #Author: Aurore Delaigle
  #Compute a version of weighted CV used in SIMEX, for binned data
  #W plays the role of the contaminated data
  #midbin is the vector of the centers of the binned data that play the role of the non contaminated data
  
  #Default values of phiU(t)=characteristic function of the errors
  #If you want to consider another error type, simply replace phiU by the characteristic function of your error type
  W=as.vector(W)
  
  if (errortype=="Lap")
  {phiU<-function(t){return(1/(1+sigU^2*t^2))}}
  if (errortype=="norm")
  {phiU<-function(t) {return(exp(-sigU^2*t^2/2))}}


  #phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure
  #you change the range of t-values, which should correspond to the support of phiK  	
  #The phiK used in this code must be the same as in the other codes (hSIMEX + NW codes)
  phiK = function(t) {return((1-t^2)^3);}
  dt = .001;
  t = seq(-1,1,dt);
  th=t/h;
  longt=length(t);
  
  
  
  
  #Compute the empirical characteristic function of W (times n) at t/h: \hat\phi_W(t/h)
  OO=matrix(rep(th,n),nrow=n,byrow=T);
  OO=matrix(rep(W,longt),ncol=longt,byrow=F)*OO;
  csO=cos(OO);
  snO=sin(OO);
  rm (OO);
  
  rehatphiW=apply(csO,2,sum);
  imhatphiW=apply(snO,2,sum);
  
  #Compute \sum_j Y_j e^{itW_j/h}
  dim(Y)=c(1,n);
  renum=Y%*%csO;
  imnum=Y%*%snO;
  
  #Compute \hat m(M_i) where M_i is the middle of the bin in which X_i (the non contaminated data) lies
  xt=matrix(rep(th,nbin),ncol=nbin,byrow=F); 
  xt=xt*matrix(rep(midbin,longt),nrow=longt,byrow=T);
  cxt=cos(xt);
  sxt=sin(xt);
  rm(xt);
  cxt=cxt[,indbin];
  sxt=sxt[,indbin];
  
  
  phiUth=phiU(th);
  matphiKU=phiK(t)/phiUth;
  dim(matphiKU)=c(1,longt)
  Den=(rehatphiW*matphiKU)%*%cxt+(imhatphiW*matphiKU)%*%sxt;
  Num=(renum*matphiKU)%*%cxt+(imnum*matphiKU)%*%sxt;
  
  
  #Compute from there the leave-one-out version \hat m_{-i}(M_i)
  
  csO=t(csO);
  snO=t(snO);
  
  Den=Den-matphiKU%*%(csO*cxt)-matphiKU%*%(snO*sxt);
  for (i in 1:n)
  {csO[,i]=csO[,i]*Y[i];
   snO[,i]=snO[,i]*Y[i];
  }
  
  Num=Num-matphiKU%*%(csO*cxt)-matphiKU%*%(snO*sxt);
  
  
  #Finally compute weighted CV, where the ith term of the sum is weighted by f_W(W_i)
  
  rhogrid=rhogrid*(2*pi*h*n)/dt;
  hW=1.06*sqrt(stats::var(W))*n^(-1/5);
  xout=outerop(W,W,"-");
  
  fWEF=rep(0,length(W))
  for  (i in 1:(length(W)))
  {for (j in 1:(length(W)))
  {xout[i,j]=stats::dnorm(xout[i,j],0,hW)}
  fWEF[i]=mean(stats::dnorm(xout[i,],0,hW));}
  
  CV=0*rhogrid;
  for (krho in 1:length(rhogrid))
  {rho=rhogrid[krho];
   dd=Den;
   dd[which(dd<rho)]=rho;
   
   mhatstar=Num/dd;
   dim(mhatstar)=c(1,n);		
   CV[krho]=sum(t(fWEF)*(Y-mhatstar)^2);
  }
  return (CV)}
