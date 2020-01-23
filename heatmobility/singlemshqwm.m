function [UF]=singlemshqwm(U0,gamma,theta,m,MU,dt,N,h,A)

l=0:(N-1);

E2=@(U) dt*laplacian5wM1d(log(U),N,h,MU)-A*dt*laplacian51d(U,N,h);
Uall=zeros(m,N);
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:));
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    S=sum(gamma(step,:));
    
    if A==0
        UF=RHS/S;
        
    else
        invmatrix=S-dt*A/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
        RHSbar=fft2(RHS);
        Ubar=RHSbar./invmatrix;
        UF=real(ifft2(Ubar));
    end
    %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    Uall(step,:)=UF;
end


