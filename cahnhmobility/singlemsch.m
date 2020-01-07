function [UF]=singlemsch(U0,gamma,theta,m,MU,dt,A,N,h)

[l,k]=meshgrid(0:(N-1));
E2=@(U) -dt*laplacian9wM(laplacian9(U,N,h)-(2*U-6*U.^2+4*U.^3),N,h,MU)+A*dt*laplacian9(laplacian9(U,N,h),N,h);

Uall=zeros(m,N,N);
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:,:));
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    S=sum(gamma(step,:));
    

    invmatrix1=sqrt(S)+1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
    invmatrix2=sqrt(S)-1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
    RHSbar=fft2(RHS);
    Ubar=(RHSbar./(invmatrix1))./(invmatrix2);
    UF=real(ifft2(Ubar));
    %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    Uall(step,:,:)=UF;
end


