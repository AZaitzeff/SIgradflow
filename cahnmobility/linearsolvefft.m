function [UN1]=linearsolvefft(U0,MU,f,S,dt,eps,N,h,C)
tol=1e-10;
maxiter=3000;
[l,k]=meshgrid(0:(N-1));

invmatrix1=sqrt(S)+1j*sqrt(C*dt*eps)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix2=sqrt(S)-1j*sqrt(C*dt*eps)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
UN=U0;

for iter=1:maxiter
    RHS=f-dt*laplacian9wM(eps*laplacian9(UN,N,h)+UN,N,h,MU)+dt*C*eps*laplacian9(laplacian9(UN,N,h),N,h);
    RHSbar=fft2(RHS);
    Ubar=(RHSbar./(invmatrix1))./(invmatrix2);
    UN1=real(ifft2(Ubar));
    if max(abs(UN1-UN))<tol
        break
    end
    UN=UN1;
end
iter