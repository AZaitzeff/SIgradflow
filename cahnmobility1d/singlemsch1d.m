function [UF]=singlemsch1d(U0,gamma,theta,m,MU,dt,eps,N,h,Force)
A2=operatormatrixpart1d(MU,N,dt,eps,h);
E2=@(U) dt*laplacian5wM1d((U.^2-1).*U+Force,N,h,MU);
Uall=zeros(m,N);
UF=U0;
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=Uall(zm-1,:);
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    S=sum(gamma(step,:));

    A=S*speye(N)+A2;
    UF=A\(RHS');

    UF=UF';
    Uall(step,:)=UF;
end