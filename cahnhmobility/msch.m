%nt=2^24;
N=2048;
mob=@(u) u.*(1-u)+1/16;
sqrtmob=@(u) sqrt(u.*(1-u)+1/16);
A=5/32;
T=20*8;
eps=1;
order=3;
nts=[2^11,2^12,2^13,2^14];
for nt=nts

[U,X,Y,h]=initializebigoval(N,eps);
Uinit=U;
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
        [U]=singlemsch(U,[[1]],[[1]],1,MU,dt,A,N,h);
    elseif order==2
        MU=mob(U);
        [Us]=singlemsch(U,[[1]],[[1]],1,MU,dt/2,A,N,h);
        MU=mob(Us);
        [U]=singlemsch(U,gamma,theta,m,MU,dt,A,N,h);
    else
        MU=mob(U);
        [Us1]=singlemsch(U,[[1]],[[1]],1,MU,dt/3,A,N,h);
        MU=mob(Us1);
        [Us2]=singlemsch(U,g2.gamma,g2.theta,g2.m,MU,dt*2/3,A,N,h);
        MU=mob(U)*1/4+mob(Us2)*3/4;
        [U]=singlemsch(U,gamma,theta,m,MU,dt,A,N,h);
    end

end
toc;
save(['../results/mschwm' num2str(order) 'order' num2str(nt) 'nt' num2str(N) 'N'],'T','dt','N','U','eps')
end