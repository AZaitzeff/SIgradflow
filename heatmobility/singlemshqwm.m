function [UF]=singlemshqwm(U0,gamma,theta,m,MU,dt,N,h,A)


E2=@(U) dt*laplacian3wM1d(log(U)-A*U,N,h,MU);
Uall=zeros(m,N);
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:));
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    
    S=sum(gamma(step,:));
    
    UF=onednonconstantlinear(MU,h,S,A*dt,RHS,N);
    %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    Uall(step,:)=UF;
end


