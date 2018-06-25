%17 April 2018
%
%Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015), 
%Robust econometric inference for stock return predictability,
%The Review of Financial Studies, 28(5), 1506-1553.
%
% The function estimates the predictive regression
% y{t}=mu+A*x{t-1}+e{t} (1), 
% and the VAR model
% x{t}=R*x{t-1}+u{t}, with R being a diagonal matrix.
% 
%
%INPUT:
% yt: a Tx1 vector which is the precticted variable y{t}

% xt: a Txkreg matrix including the regressors as columns. 
% The program automatically includes an intercept in the predictive
% regression, so there is no need for a column of 1's.
% There is no need to input the lagged series; the program lags the series.
%
% K: is the horizon (special case K=1 corresponds to a short-horizon
% regression)
%
% printres: set printres=1 if you want the results to be printed;

%OUTPUT:
% Aols: a (kreg+1)x1 vector which is the OLS estimator of the intercept mu
% (1st element) and the slope coefficients (A) 
%  
% Aivx: a 1xkreg vector which is the IVX estimator of A in (1)
%
% Wivx : is a 2x1 vector with the first element being the Wald statistic of 
% a test for overall sifnificance, 
% whereas the second gives the corresponding p-value (from the
% chi-square distribution)
%
% WivxInd: is a 2xkreg matrix 
% the first row gives the individual test of siginificance for each
% predictor (i.e. restricting each predictor's coefficient to be equal to
% 0 letting the rest coefficients free to be estimated), 
% whereas the second row gives the corresponding p-values (for a two-sided test).
%
% Q is the variance-covariance matrix of the IVX estimator for A.
%
% corrmat gives the correlation between the residuals of the predictive
% regression (e{t}) and the residuals of the autoregressive part (u{t})


function [Aols,Aivx,Wivx,WivxInd,Q,corrmat]=ivxlh(yt,xt,K,printres)
xlag=xt(1:end-1,:);
xt=xt(2:end,:);
y=yt(2:end,:);

[nn,l]=size(xlag); 
X=[ones(nn,1) xlag];

Wivx=zeros(2,1);
WivxInd=zeros(2,l);

%predictive regression residual estimation 
[Aols,bhat,epshat]=regress(y,X); 

rn=zeros(l,l);
for i=1:l
    rn(i,i)=regress(xt(:,i),xlag(:,i));
end

%autoregressive residual estimation 
u=xt-xlag*rn;

%residuals' correlation matrix
corrmat=corrcoef([epshat u]);

%covariance matrix estimation (predictive regression)
covepshat=epshat'*epshat/nn;
covu=zeros(l,l);
for t=1:nn
    covu=covu+u(t,:)'*u(t,:);
end

%covariance matrix estimation (autoregression)
covu=covu/nn;
covuhat=zeros(1,l);
for i=1:l
    covuhat(1,i)=sum(epshat.*u(:,i));
end

%covariance matrix between 'epshat' and 'u'
covuhat=covuhat'/nn; 

m=floor(nn^(1/3)); 
uu=zeros(l,l);
for h=1:m
    a=zeros(l,l);
    for t=(h+1):nn
        a=a+u(t,:)'*u(t-h,:);
    end
    uu=uu+(1-h/(m+1))*a;
end 
uu=uu/nn;
Omegauu=covu+uu+uu'; 

q=zeros(m,l);
for h=1:m
    p=zeros(nn-h,l);
    for t=(h+1):nn
        p(t-h,:)=u(t,:)*epshat(t-h)';
    end
    q(h,:)=(1-h/(1+m))*sum(p);
end
residue=sum(q)/nn;
Omegaeu=covuhat+residue'; 

%instrument construction
n=nn-K+1;
Rz=(1-1/(nn^0.95))*eye(l); 
diffx=xt-xlag; 
z=zeros(nn,l);
z(1,:)=diffx(1,:);
for i=2:nn
    z(i,:)=z(i-1,:)*Rz+diffx(i,:);
end
Z=[zeros(1,l);z(1:n-1,:)];


zz=[zeros(1,l);z(1:nn-1,:)];
ZK=zeros(n,l);
for i=1:n
    ZK(i,:)=sum(zz(i:i+K-1,:),1);
end

yy=zeros(n,1);
for i=1:n
    yy(i)=sum(y(i:i+K-1));
end 
xK=zeros(n,l);
for i=1:n
    xK(i,:)=sum(xlag(i:i+K-1,:),1);
end 

meanxK=mean(xK);
Yt=yy-mean(yy);
Xt=zeros(n,l);
for i=1:l
    Xt(:,i)=xK(:,i)-meanxK(:,i)*ones(n,1);
end

Aivx=Yt'*Z*pinv(Xt'*Z);
meanzK=mean(ZK);

FM=covepshat-Omegaeu'*Omegauu^(-1)*Omegaeu;
M=ZK'*ZK*covepshat-n*meanzK'*meanzK*FM;

H=eye(l);
Q=H*pinv(Z'*Xt)*M*pinv(Xt'*Z)*H';
Wivx(1,1)=(H*Aivx')'*pinv(Q)*(H*Aivx');   
Wivx(2,1)= 1-cdf('chi2',Wivx(1,1), l);
    
WivxInd(1,:)=Aivx./((diag(Q)).^(1/2))';
WivxInd(2,:)=1-cdf('chi2',WivxInd(1,:).^2, 1);
    
    if printres==1;
    disp('Horizon: ');
    disp(K);
    disp('               ');
    disp('IVX estimator: ');
    disp(Aivx);
    disp('               ');
    disp('Individual Tests of significance (1st row is the ): ');
    disp(WivxInd);
    disp('               ');
    disp('Test of overall significance: ');
    disp(Wivx);
    disp('Residuals correlation matrix: ');
    disp(rescorrmat);
    end;


end

