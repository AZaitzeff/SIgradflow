nt=2^24;
N=1024;

%N=512;
%nt=2^15;

%N=256;
%nt=2^10;

T=1/2;
eps=0.0025;
theta=pi/4;
%[U,X,Y,h]=initializecos(N);
[U,X,Y,h]=initializebigoval(N,eps,theta);

Uinit=U;
dt=T/nt;
mob=@(u) (1-sqrt(eps))*(1-u.^2).^2+sqrt(eps);
sqrtmob=@(u) sqrt((1-u.^2));
Wp=@(u) (u.^2-1).*u;

%Wp=@(u) 1/4*log((1-u)./(1-u))-u;
A=1;
[l,k]=meshgrid(0:(N-1));
invmatrix01=1+1j*sqrt(eps*dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix02=1-1j*sqrt(eps*dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
%invmatrix1=sqrt(3/2)+1j*sqrt(eps*dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
%invmatrix2=sqrt(3/2)-1j*sqrt(eps*dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));


%invmatrix01=1+1j*sqrt(eps*dt*A)/(h^2)*(-4+2*cos(2*pi*l/N)+2*cos(2*pi*k/N));
%invmatrix02=1-1j*sqrt(eps*dt*A)/(h^2)*(-4+2*cos(2*pi*l/N)+2*cos(2*pi*k/N));
%invmatrix1=sqrt(3/2)+1j*sqrt(eps*dt*A)/(h^2)*(-4+2*cos(2*pi*l/N)+2*cos(2*pi*k/N));
%invmatrix2=sqrt(3/2)-1j*sqrt(eps*dt*A)/(h^2)*(-4+2*cos(2*pi*l/N)+2*cos(2*pi*k/N));

tic;
Upre=U;
for t=1:nt
    
    %if t==1
        
        F=U-dt*laplacian9wM(laplacian9(eps*U,N,h)-Wp(U),N,h,mob(U))+A*eps*dt*laplacian9(laplacian9(U,N,h),N,h);
        Fbar=fft2(F);
        Ubar=(Fbar./(invmatrix01))./(invmatrix02);
        U=real(ifft2(Ubar));
%     else
%         
%         F=2*U-1/2*Upre+...
%             2*(-dt*laplacian5wM(laplacian5(eps*U,N,h)-Wp(U),N,h,mob(U))+eps*A*dt*laplacian5(laplacian5(U,N,h),N,h))-...
%         (-dt*laplacian5wM(laplacian5(eps*Upre,N,h)-Wp(Upre),N,h,mob(Upre))+eps*A*dt*laplacian5(laplacian5(Upre,N,h),N,h));
%         Upre=U;
%         Fbar=fft2(F);
%         Ubar=(Fbar./(invmatrix1))./(invmatrix2);
%         U=real(ifft2(Ubar));
%         
%     end
    if mod(t,2^15)==0 && t~=nt
        save(['../results/siechwm' num2str(t) 't' num2str(nt) 'N' num2str(N)],'T','dt','N','U','Upre','eps','X','Y','h')
    end
end
toc;
save(['../results/siechwm' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps','X','Y','h')