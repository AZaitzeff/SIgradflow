
N=4096+1;
T=1;
N=512;
%nt=2^18;


nts=[2^18];

numofts=size(nts,2);

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
        Ubar=dt*laplacian51d(nthroot(U,3).^(5),N,h)+U;
        U=dt/2*laplacian51d(nthroot(Ubar,3).^(5),N,h)+1/2*(U+Ubar);
        %der(t)=max(abs(firstder51d(U',N,h)));
    end
    toc;
    save(['2pc' num2str(N) 'N' num2str(nt) 'nt.mat'],'U','x','h')
end