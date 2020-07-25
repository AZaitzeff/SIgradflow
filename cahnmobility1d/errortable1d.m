
vars=load(['../results/siechwm1d4194304N2048.mat']);
%vars=load(['../results/siechwm1d1048576N2048.mat']);
trueU=vars.U;
N=2048;
eps=.05^2;
[~,trueX,h]=initialize1d(N,eps);
%trueU=vars.U(1:end,1:end);

nts=[2^7,2^8,2^9,2^10,2^11];
%nts=[2^8,2^9,2^10,2^11,2^12,2^13];
n=numel(nts);
error1=zeros(2,n);
Ns=[2048,2048,2048,2048,2048,2048,2048];
order=3;
%load(['../results/mschwm' num2str(order) 'order' num2str(2^6) 'nt' num2str(N) 'N']);
for numt=1:n
    %Up=U;
    nt=nts(numt);
    N=Ns(numt);
    %N=ceil(512*2^((numt-1)/4));
    %N=128;
    [~,X,h]=initialize1d(N,eps);
    load(['../results/mschwm1d' num2str(order) 'order' num2str(nt) 'nt' num2str(N) 'N']);
    tU= interp1(trueX,trueU,X,'spline');
    %figure
    %imagesc(U)
    error1(1,numt)=sqrt(h*sum((U(:)-tU(:)).^2));
    %error1(1,numt)=max(abs(U(:)-tU(:)));
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end