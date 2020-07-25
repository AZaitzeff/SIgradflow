%nt=2^24;
%N=2048*4;
%nt=2^28;

N=2048;
nt=2^22;

%N=128;
%nt=2^10;

T=1/8;
eps=0.0025;%really epsilon squared
[U,X,h]=initialize1d(N,eps);
Uinit=U;
dt=T/nt;
Force=tanh(1/(10*sqrt(eps))*cos(2*pi*X));
%mob=@(u) (1-u.^2).^2;
mob=@(u) (1-sqrt(eps))*(1-u.^2).^2+sqrt(eps);
sqrtmob=@(u) (1-u.^2);
Wp=@(u) (u.^2-1).*u;

%Wp=@(u) 1/4*log((1-u)./(1-u))-u;
A=1;
l=0:(N-1);
invmatrix01=1+1j*sqrt(eps*dt*A)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
invmatrix02=1-1j*sqrt(eps*dt*A)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
invmatrix1=sqrt(3/2)+1j*sqrt(eps*dt*A)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
invmatrix2=sqrt(3/2)-1j*sqrt(eps*dt*A)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));



tic;
Upre=U;
for t=1:nt
    if t==1
        F=U-dt*laplacian5wM1d(laplacian51d(eps*U,N,h)-Wp(U)-Force,N,h,mob(U))+A*eps*dt*laplacian51d(laplacian51d(U,N,h),N,h);
        Fbar=fft(F);
        Ubar=(Fbar./(invmatrix01))./(invmatrix02);
        U=real(ifft(Ubar));
    else

        F=2*U-1/2*Upre+...
            2*(-dt*laplacian5wM1d(laplacian51d(eps*U,N,h)-Wp(U)-Force,N,h,mob(U))+A*eps*dt*laplacian51d(laplacian51d(U,N,h),N,h))-...
            (-dt*laplacian5wM1d(laplacian51d(eps*Upre,N,h)-Wp(Upre)-Force,N,h,mob(Upre))+A*eps*dt*laplacian51d(laplacian51d(Upre,N,h),N,h));
        Upre=U;
        Fbar=fft(F);
        Ubar=(Fbar./(invmatrix1))./(invmatrix2);
        U=real(ifft(Ubar));
    end
end
toc;
%tU= interp1(trueX,trueU,x,'spline');
%sqrt(h*sum((U(:)-tU(:)).^2))
save(['../results/siechwm1d' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps','X','h')