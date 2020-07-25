addpath('..')
N=4096*4+1;
T=1;
order=3;
%nts=[2^21];
nts=[2^4,2^5,2^6,2^7,2^8];%,2^8,2^9];
numofts=size(nts,2);


error1=zeros(2,numofts);
var=load('trbdf216385N65536nt.mat');
trueU=var.U;
[U,x,h] = initializegauss(N,2/3);
%U = conv(U,guass,'same')*h;

Uinit=U;
mob=@(U) U;
%eps=[4e-5,;
energy=zeros(1,2^5+1);
energy(1)=sum(3/2*U.^(5/3));
for iter=1:numofts
nt=nts(iter);
U=Uinit;
dt=T/nt;




if order==2
    [gamma,m]=getgamma(2,0);
    
else
    [si1125gamma,si1125m]=getgamma(2,1);
    [si3ordergamma,si3orderm]=getgamma(3,0);
    
end
tic;
for t=1:nt
    if order==1
        MU=mob(U);
        [U]=singlemspmewm53(U,[[1]],1,MU,dt,N,h);
        energy(t+1)=sum(3/2*U.^(5/3));
    elseif order==2
        MU=mob(U);
        [Us]=singlemspmewm53(U,[[1]],1,MU,dt/2,N,h);
        MU=mob(Us);
        [U]=singlemspmewm53(U,gamma,m,MU,dt,N,h);
        energy(t+1)=sum(3/2*U.^(5/3));
    elseif order==3
        MUn=mob(U);
        
        [U1h]=singlemspmewm53(U,[[1]],1,MUn,dt/6,N,h);
        MU1h=mob(U1h);
        
        [U2h1s]=singlemspmewm53(U,[[1]],1,MUn,2*dt/5,N,h);
        MU2h1s=mob(U2h1s);
        [U2h]=singlemspmewm53(U,si1125gamma,si1125m,MU2h1s,dt*5/6,N,h);
        MU2h=mob(U2h);
        
        [U1]=singlemspmewm53(U,si3ordergamma,si3orderm,MU1h,dt/2,N,h);
        [U]=singlemspmewm53(U1,si3ordergamma,si3orderm,MU2h,dt/2,N,h);
        energy(t+1)=sum(3/2*U.^(5/3));
    end
    
end
toc;
%U=U-eps;
%save(['pmewm' num2str(order) 'o' num2str(N) 'N' num2str(nt) 'nt.mat'],'U','x','h','energy','T','dt')
%error1(1,iter)=sqrt(var.h*sum((U(1:step:end)-trueU).^2));
%plot(x,U);
error1(1,iter)=sqrt(h*sum((U-trueU).^2));
end
% 
for i=2:numofts
      error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end