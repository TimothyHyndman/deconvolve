function hPI = PI_deconvUesth3(W1,W2)

%W1 and W2: replicated data from contaminated sample. W1 and W2 are of the same size

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%								WARNINGS:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% The range of t-values -1 and 1 correspond to the support of phiK. 
% If you change phiK and take a kernel for which phiK is not supported on [-1,1] you have to change -1 and 1 accordingly
%
% muK2 below is the second moment \int x^2 K(x) dx  of the kernel K defined below through phiK. 
% If you change phiK, you MUST change muK2 accordingly. 
% !!! The phiK here must match the phiK used to compute the estimator itself.!!! 
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------


phiK = @(t) (1-t.^2).^3;
deltat = .0002;
dt = .0002;
t = (-1:dt:1)';
t=reshape(t,length(t),1);
muK2 = 6;


%------------------------------------------------------
% Compute bandwidth for density that is a multiple of g
%------------------------------------------------------
n=length(W1);
W1=reshape(W1,n,1);
W2=reshape(W2,n,1);
W=[W1;W2];
n=2*n;

%grid of h values where to search for a solution
maxh=(max(W)-min(W))/10
hnaive=1.06*sqrt(var(W))*n^(-1/5);
hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;
lh = length(hgrid);

%------------------------------------------------------------------------------------------
%Estimate phiU on a large grid (the one corresponding to the smallest bandwidth considered)
%------------------------------------------------------------------------------------------

tU=t/hgrid(1);
Diff=W1-W2;
Diff=Diff((W1~=0)&(W2~=0));
nDiff=length(Diff);
tout=outerop(tU,Diff,'*');
phiU=sum(cos(tout),2)/nDiff;
phiU(phiU<0)=0;
phiU=sqrt(abs(phiU));


%Find the smallest and largest value of t for which we use this estimator, i.e. find where it becomes unreliable. For this find where  |phiU|< n^(-1/4). 
tmp=tU(abs(phiU)<nDiff^(-0.25));
t2=min(tmp(tmp>0));

%Use a Laplace approximation elsewhere and then use a spline approximation of phiU to put the pieces together.
tlim=[-t2,t2]
ppphiU=spline(tU,phiU);
hatvarU=var(Diff)/2;

%-----------------------------------------------
%Continue with the bandwidth selection procedure
%-----------------------------------------------

stdevx = max(sqrt(var(W) - hatvarU),1/n); % std(X)
th3 = 0.5289277*stdevx^(-7); % Estimate theta4 by NR     
hgrid=reshape(hgrid,1,lh);
toverh=t*(1./hgrid);

phiK2=(phiK(t)).^2;
%Use a Laplace approximation in the tails, as done in DH (2016), JRSSB and also in DHM (2008)
phiU2=phiUspline(toverh,hatvarU,tlim,ppphiU).^2;



rr=2;
% Find h2 for th2
term1= -hgrid.^2*muK2*th3;
term2=repmat(t.^(2*rr).*phiK2, 1, lh)./phiU2;
term2=sum(term2,1)*dt./(2*pi*n*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;
indh2=find(ABias2==min(ABias2),1,'first')
h2 = hgrid(indh2);


OO=outerop(t/h2,W,'*');
%Estimate empirical characersitic fucntion of W
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
clear OO;
normhatphiW2=rehatphiW.^2+imhatphiW.^2;
th2 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh2))*dt/(2*pi*h2^(2*rr+1));


term1=hgrid.^4*muK2^2*th2/4;
term2=repmat(phiK2,1,lh)./phiU2;
term2=sum(term2,1)*dt./(2*pi*n*hgrid);
AMISE=term1+term2;

indh=find(AMISE==min(AMISE),1,'first')
hPI = hgrid(indh);

