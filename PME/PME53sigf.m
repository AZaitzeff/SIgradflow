addpath('..')
if order==2
    load('si2nd.mat')
elseif order==3
    load('si3rd.mat')
end
fac=3/4;
N=4096*4+1;
T=1;

nts=2.^(11:17);
numofts=size(nts,2);
MU=ones(N,1);
error1=zeros(2,numofts);
var=load('trbdf216385N65536nt.mat');
trueU=var.U;
[U,x,h] = initializegauss(N,2/3);
Uinit=U;
for numt=1:numofts

    %N=ceil(512*2^((numt-2)/2));
    Uall=zeros(m+1,N);
    nt=nts(numt);
    U=Uinit;
    dt=T/nt;
    tic;
for t=1:nt
    curU=U;
    Uall(1,:)=U;
    for step=1:m
        Utemp=squeeze(Uall(1,:));
        Un=Utemp*gamma(step,1)+dt*theta(step,1)*laplacian31d(nthroot(Utemp,3).^(5)-fac*Utemp,N,h);

        for zm=2:step
            Utemp=squeeze(Uall(zm,:));
            Un=Un+Utemp*gamma(step,zm)+dt*theta(step,zm)*laplacian31d(nthroot(Utemp,3).^(5)-fac*Utemp,N,h);
        end

        S=sum(gamma(step,:));
        
        [curU,~]=onedsolvematrix(MU,MU',h,S,fac*dt,Un,N);
        
        Uall(step+1,:)=curU;

    end

    U=curU;

    

end

toc;

error1(1,numt)=sqrt(h*sum((U-trueU).^2));

end

for i=2:numofts
      error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end