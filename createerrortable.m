%nt=1024*2^6;
vars=load(['results/si2ac1048576N2048']);
trueU=vars.U;
N=2048;
[~,trueX,trueY,h]=initializebigcircle(N,1);
%[~,r,~]=initializebigcirclepolar(4096,eps);
nts=2.^(12:16);
n=numel(nts);
error1=zeros(2,n);
for numt=1:n
    nt=nts(numt);
    N=2048;
    [~,X,Y,h]=initializebigcircle(N,1);
    load(['results/multistepac2s' num2str(nt) 'N' num2str(N)]);
    tU= interp2(trueX,trueY,trueU,X,Y,'spline');

    error1(1,numt)=sqrt(sum((tU(:)-U(:)).^2)*h^2);
end
    %plot(trueX(512,:),trueZ(512,:),'k');hold on;
for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end