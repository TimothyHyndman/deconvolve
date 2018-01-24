function fXdecUK=gdecUunknown(xx,W1,W2,h)

%Author: Aurore Delaigle
%Compute the deconvolution kernel density estimator

%xx: vector of x-values where to compute the deconvolution kernel density estimator
%W1 and W2: replicated data from contaminated sample. W1 and W2 are of the same size
%h bandwidth

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%								WARNINGS:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% The range of t-values -1 and 1 correspond to the support of phiK. 
% If you change phiK and take a kernel for which phiK is not supported on [-1,1] you have to change -1 and 1 accordingly
% !!! The phiK here must match the phiK used to compute the bandwidth (PI, CV or other).!!! 
%
% The DKDE can also be computed using the Fast Fourier Transform, which is a bit more complex. 
% See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution.   Statistics and Computing,  17,  349 - 355
% However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------



%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
%you change the range of t-values, which should correspond to the support of phiK
phiK = @(t) (1-t.^2).^3;


%Range of t-values (must correspond to the domain of phiK)
%If you get poor results, check that this grid is fine enough
deltat = .0002;
tt = (-1:deltat:1)';
tt=reshape(tt,length(tt),1);

%----------------------------------------
%Estimate phiU using the replicated data
%----------------------------------------

Diff=W1-W2;
nDiff=length(Diff);
tout=outerop(tt/h,Diff,'*');
phiU=sum(cos(tout),2)/nDiff;
phiU(phiU<0)=0;
phiU=sqrt(abs(phiU));

tU=tt/h;
%Find the smallest and largest value of t for which we use this estimator, i.e. find where it becomes unreliable. For this find where  |phiU|< n^(-1/4). 
tmp=tU(abs(phiU)<nDiff^(-0.25));
t2=min(tmp(tmp>0));
if(length(t2)==0)
	t2=tU(length(tU));
end

%Use a Laplace approximation elsewhere and then use a spline approximation of phiU to put the pieces together.
tlim=[-t2,t2];
ppphiU=spline(tU,phiU);
hatvarU=var(Diff)/2;
phiUS=phiUspline(tt/h,hatvarU,tlim,ppphiU);

%--------------------
% Estimate FT of fX. 
%--------------------

W1=reshape(W1,nDiff,1);
W2=reshape(W2,nDiff,1);
WTout=[W1;W2];
ntout=2*nDiff;

OO=outerop(tt/h,WTout,'*');
rehatphiW=sum(cos(OO),2)/ntout;
imhatphiW=sum(sin(OO),2)/ntout;
rehatphig=rehatphiW./phiUS;
imhatphig=imhatphiW./phiUS;

xt=outerop(tt/h,xx,'*');
longx=length(xx);
fXdecUK=cos(xt).*repmat(rehatphig,1,longx)+sin(xt).*repmat(imhatphig,1,longx);
fXdecUK=sum(fXdecUK.*repmat(phiK(tt),1,longx),1)/(2*pi)*deltat/h;
fXdecUK(fXdecUK<0)=0*fXdecUK(fXdecUK<0);

