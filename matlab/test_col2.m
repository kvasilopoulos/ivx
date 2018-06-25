data = xlsread('monthly.xlsx',1);
K = 1;
yt = data(:, 13); 
xt = data(:, 2);

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