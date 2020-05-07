function [UF,RHS]=singlemsch(U0,gamma,theta,m,MU,dt,eps,N,h)
%A2=operatormatrixpart(MU,N,dt,eps,h);
tol1=1e-6;
E2=@(U) dt*laplacian5wM((U.^2-1).*U,N,h,MU);
Uall=zeros(m,N,N);
UF=U0;
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:,:));
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    S=sum(gamma(step,:));
    [UF]=iterspecialsolve(UF,MU,dt,eps,N,h,S,RHS,tol1);
    Uall(step,:,:)=UF;
end


