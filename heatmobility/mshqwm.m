
N=8192*2+1;
T=1/10;
order=3;
nts=[2^3,2^4,2^5,2^6,2^7,2^8];
%nts=[2^11,2^12,2^13];
A=1;
numofts=size(nts,2);
error1=zeros(2,numofts);

[U,x,h]=initializesin(N);
Uinit=U;
mob=@(U) U;
trueU=cos(pi*x)*exp(-T*(pi^2))+2;
for iter=1:numofts
nt=nts(iter);
U=Uinit;
dt=T/nt;




if order==2
    load('../si2nd.mat')
    
else
    si1=load('../scheme1si.mat');
    si1125=load('../scheme1125si.mat');
    si3order=load('../si3rd.mat');
    
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
        MUn=mob(U);
        
        [U1h]=singlemshqwm(U,si1.gamma,si1.theta,si1.m,MUn,dt/6,N,h,A);
        MU1h=mob(U1h);
        
        [U2h1s]=singlemshqwm(U,[[1]],[[1]],1,MUn,2*dt/5,N,h,A);
        MU2h1s=mob(U2h1s);
        [U2h]=singlemshqwm(U,si1125.gamma,si1125.theta,si1125.m,MU2h1s,dt*5/6,N,h,A);
        MU2h=mob(U2h);
        
        [U1]=singlemshqwm(U,si3order.gamma,si3order.theta,si3order.m,MU1h,dt/2,N,h,A);
        [U]=singlemshqwm(U1,si3order.gamma,si3order.theta,si3order.m,MU2h,dt/2,N,h,A);
    end
    
end
toc;
error1(1,iter)=sqrt(h*sum((U-trueU).^2));
end

for i=2:numofts
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end