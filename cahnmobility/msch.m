%nt=2^24;
Ns=[256,ceil(256*sqrt(2)),ceil(256*2)];
mob=@(u) (1-u.^2);
A=0;
T=1/10;
eps=3.2e-4;
order=2;
nts=[2^7,2^8,2^9,2^10,2^11];

for nt=nts
    for N=Ns

[~,X,Y,h]=initializerandom(N);
init=load('init256.mat');
U=interp2(init.X,init.Y,init.U,X,Y,'spline');
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
        [U]=singlemsch(U,1,1,MU,dt,eps,N,h,A);
        %sum(U.^2)
    elseif order==2
        MU=mob(U);
        [Us]=singlemsch(U,1,1,MU,dt/2,eps,N,h,A);
        MU=mob(Us);
        [U]=singlemsch(U,gamma,m,MU,dt,eps,N,h,A);
        %sum(U.^2)
    elseif order==3
        MUn=mob(U);
        
        [U1h]=singlemsch(U,1,1,MUn,dt/6,eps,N,h,A);
        MU1h=mob(U1h);
        
        [U2h1s]=singlemsch(U,1,1,MUn,2*dt/5,eps,N,h,A);
        MU2h1s=mob(U2h1s);
        [U2h]=singlemsch(U,si1125gamma,si1125m,MU2h1s,dt*5/6,eps,N,h,A);
        MU2h=mob(U2h);
        
        [U1]=singlemsch(U,si3ordergamma,si3orderm,MU1h,dt/2,eps,N,h,A);
        [U]=singlemsch(U1,si3ordergamma,si3orderm,MU2h,dt/2,eps,N,h,A);
        %sum(U.^2)
    end
    
end
toc;
save(['../results/mschwm' num2str(order) 'order' num2str(nt) 'nt' num2str(N) 'N'],'T','dt','N','U','eps')
end
end