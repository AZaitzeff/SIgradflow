%nt=2^24;
%N=2048;
N=512;
T=1/50;
order=3;
nts=[2^6,2^7,2^8,2^9,2^10];
%nts=[2^11,2^12,2^13];
A=1;
numofts=size(nts,2);
error1=zeros(2,numofts);

[U,x,h]=initializesin(N);
Uinit=U;
mob=@(U) U;
trueU=sin(2*pi*x)*exp(-T*(4*pi^2))+2;
for iter=1:numofts
nt=nts(iter);
[u,x,h]=initializesin(N);
U=Uinit;
dt=T/nt;




if order==2
    load('../si2nd.mat')
    
else
    g2=load('../si2nd.mat');
    load('../si3rd.mat')
    
end

tic;
for t=1:nt
    if order==1
        MU=mob(U);
        [U]=singlemshqwm(U,[[1]],[[1]],1,MU,dt,N,h,A);
    elseif order==2
        MU=mob(U);
        [Us]=singlemshqwm(U,[[1]],[[1]],1,MU,dt/2,N,h,A);
        MU=mob(Us);
        [U]=singlemshqwm(U,gamma,theta,m,MU,dt,N,h,A);
    elseif order==3
        MU1=mob(U);
        [Us1]=singlemshqwm(U,[[1]],[[1]],1,MU1,dt/3,N,h,A);
        MU2=mob(Us1);
        [Us2]=singlemshqwm(U,g2.gamma,g2.theta,g2.m,MU2,dt*2/3,N,h,A);
        MU3=mob(U)*1/4+mob(Us2)*3/4;
        [U]=singlemshqwm(U,gamma,theta,m,MU3,dt,N,h,A);
    end
    
end
toc;
error1(1,iter)=sqrt(sum((U-trueU).^2));
end

for i=2:numofts
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end