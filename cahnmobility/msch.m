%nt=2^24;
Ns=[1024,1024,1024,1024,1024,1024];
eps=0.0025;
mob=@(u) (1-sqrt(eps))*(1-u.^2).^2+sqrt(eps);
mobpp=@(u) -4*(1-sqrt(eps))*(1-3*u.^2);
Wp=@(u) (u.^2-1).*u;
T=1/2;


angle=pi/4;
order=3;
nts=[2^4,2^5,2^6,2^7,2^8,2^9];
n=numel(nts);
for i=1:n

nt=nts(i);
N=Ns(i);

[U,X,Y,h]=initializebigoval(N,eps,angle);
%[U,X,Y,h]=initializecos(N);
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
        [U]=singlemsch(U,[[1]],[[1]],1,MU,dt,eps,N,h);
    elseif order==2
        MU=mob(U);
        [Us]=singlemsch(U,[[1]],[[1]],1,MU,dt/2,eps,N,h);
        MU=mob(Us);
        [U]=singlemsch(U,gamma,theta,m,MU,dt,eps,N,h);
    elseif order==3
        MUn=mob(U);
        
        [U1h]=singlemsch(U,si1.gamma,si1.theta,si1.m,MUn,dt/6,eps,N,h);
        MU1h=mob(U1h);
        w1=laplacian5wM(laplacian5(eps*U1h,N,h)-Wp(U1h),N,h,MU1h);
        fullMU1h=MU1h-1/72*k^2*mobpp(U1h).*w1.^2;
        
        [U2h1s]=singlemsch(U,[[1]],[[1]],1,MUn,2*dt/5,eps,N,h);
        MU2h1s=mob(U2h1s);
        
        
        [U2h]=singlemsch(U,si1125.gamma,si1125.theta,si1125.m,MU2h1s,dt*5/6,eps,N,h);
        MU2h=mob(U2h);
        w2=laplacian5wM(laplacian5(eps*U2h,N,h)-Wp(U2h),N,h,MU2h);
        fullMU2h=MU2h-1/72*k^2*mobpp(U2h).*w2.^2;
        
        
        [U1]=singlemsch(U,si3order.gamma,si3order.theta,si3order.m,fullMU1h,dt/2,eps,N,h);
        [U]=singlemsch(U1,si3order.gamma,si3order.theta,si3order.m,fullMU2h,dt/2,eps,N,h);
    end
    
end
toc;
save(['../results/mschwm' num2str(order) 'order' num2str(nt) 'nt' num2str(N) 'N'],'T','dt','N','U','eps')
end