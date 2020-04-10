%addpath('../findgamma/')
%[gamma,m]=getgamma(3,0);
load('si3rd.mat')
T=2;
trueval=2*atan(exp(-T));
Ns=[16,32,64,128,256,512,1024,2048,4096];
maxiter=1000;
n=numel(Ns);
error1=zeros(2,n);
for iter=1:n
    N=Ns(iter);
    k=T/N;
    t=0:k:T;
    u=zeros(N,m+1);
    
    for i=1:N
        if i==1
          u(1,1)=pi/2;  
        else
          u(i,1)=u(i-1,m+1);
        end
        %u(i+1)=u(i)/(1+k);
        for j=2:m+1
            ubar=sum(gamma(j-1,1:j-1).*u(i,1:j-1));
            S=sum(gamma(j-1,1:j-1));
            curu=1/(S+2*k)*(-k*sum((sin(u(i,1:j-1))-2*u(i,1:j-1)).*theta(j-1,1:j-1))+ubar);
            u(i,j)=curu;
        end
    end

    error1(1,iter)=abs(u(N,m+1)-trueval);
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end