function [UF]=singlemspmewm(U0,gamma,m,MU,dt,N,h)

Uall=zeros(m,N);

for step=1:m
    RHS=U0*gamma(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:));
        RHS=RHS+curU*gamma(step,zm);
    end
    S=sum(gamma(step,:));
    
    UF=onednonconstantlinear(MU,h,S,2*dt,RHS,N);

    Uall(step,:)=UF;
end


