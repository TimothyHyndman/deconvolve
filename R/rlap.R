#Generate a random matrix of size n1 x n2 from a Laplace(sigLap). Note that the variance of a Laplace(sigLap) is equal to 2 sigLap^2. Thus sigLap is NOT the standard deviation
rlap=function(sigLap,n1,n2){
  y=matrix(runif(n1*n2,0,1),nrow=n1,ncol=n2,byrow=T);
  reponse= sigLap*log(2*y);
  reponse[which(y>0.5)]=-sigLap*log(2-2*y[which(y>0.5)]);
  return(reponse)}
