addpath('..')
N=4096*4+1;
T=1;
%nt=2^18;
maxiter=100;
tol=1e-10;
MU=ones(N,1);
nts=[2^15,2^16];
numofts=size(nts,2);
%error1=zeros(2,numofts);
%var=load('2pc512N262144nt.mat');
%trueU=var.U;

%numofts=size(nts,2);

%[U,x,h] = initializebump(N);
[U,x,h] = initializegauss(N,2/3);
Uinit=U;


for iter=1:numofts
    nt=nts(iter);
    U=Uinit;
    dt=T/nt;
    der=zeros(1,nt);
    tic;
    for t=1:nt
        RHS=dt/4*laplacian31d(nthroot(U,3).^(5),N,h)+U;
        UN12=U;
    for i=1:maxiter
        val=RHS-(UN12-dt/4*laplacian31d(nthroot(UN12,3).^5,N,h));
        
        if sqrt(sum(abs(val).^2))<tol
            break
        end
        [delU,~]=onedsolvematrix(MU,nthroot(UN12,3).^2,h,1,5/3*dt/4,val,N);
        UN12=UN12+delU;
        
    end
    
        RHS=(4*UN12-U)/3;
        U=UN12;
    for i=1:maxiter
        val=RHS-(U-dt/3*laplacian31d(nthroot(U,3).^5,N,h));
        
        if sqrt(sum(abs(val).^2))<tol
            break
        end
        [delU,~]=onedsolvematrix(MU,nthroot(U,3).^2,h,1,5/3*dt/3,val,N);
        U=U+delU;
        
    end
    end
    toc;
    %error1(1,iter)=sqrt(h*sum((U-trueU).^2));
    save(['trbdf2' num2str(N) 'N' num2str(nt) 'nt.mat'],'U','x','h')
end

% for i=2:numofts
%       error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
% end