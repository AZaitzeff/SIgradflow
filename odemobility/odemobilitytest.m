order=3;

if order==2
    load('../si2nd.mat')
    
else
    g2=load('../si2nd.mat');
    load('../si3rd.mat')
    
end
T=1;
Ns=[128,256,512,2^10,2^11,2^12];

A=1;
%truex=exp(-5*T);
%truey=2*exp(5*T);

truex=1/(-4+5*exp(4*T))^(1/4);
truey=2/(-4+5*exp(4*T))^(1/4);
n=size(Ns,2);
error1=zeros(2,6);
for iter=1:n
    N=Ns(iter);
    k=T/N;
    t=0:k:T;
    ux=zeros(1,N+1);
    uy=zeros(1,N+1);
    
    for i=1:N
        
        if i==1
          ux0=1;
          uy0=2;
          ux(1)=ux0;
          uy(1)=uy0;
        else
          ux0=uxf;
          uy0=uyf;
        end
        if order==2
            MU=1+ux0^2*uy0^2;
            [sux,suy]=singlemultistep(ux0,uy0,[[1]],[[1]],1,MU,k/2,A);
            MU=1+sux^2*suy^2;
            [uxf,uyf]=singlemultistep(ux0,uy0,gamma,theta,m,MU,k,A);
        else
            MU=1+ux0^2*uy0^2;
            [sux,suy]=singlemultistep(ux0,uy0,[[1]],[[1]],1,MU,k/3,A);
            MU=1+sux^2*suy^2;
            [uxs,uys]=singlemultistep(ux0,uy0,g2.gamma,g2.theta,g2.m,MU,k*2/3,A);
            MU=(1+ux0^2*uy0^2)*1/4+(1+uxs^2*uys^2)*3/4;
            [uxf,uyf]=singlemultistep(ux0,uy0,gamma,theta,m,MU,k,A);
        end
        
        ux(i+1)=uxf;
        uy(i+1)=uyf;
    end

    error1(1,iter)=sqrt((uxf-truex)^2+(uyf-truey)^2);
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end