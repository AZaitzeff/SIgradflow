%nt=2^24;
N=2048;
nt=2^19;


T=20*8;
eps=1;
[U,X,Y,h]=initializebigoval(N,eps);
Uinit=U;
dt=T/nt;
mob=@(u) u.*(1-u)+1/16;
sqrtmob=@(u) sqrt(u.*(1-u)+1/16);
A=8/32;
[l,k]=meshgrid(0:(N-1));
invmatrix01=1+1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix02=1-1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix1=sqrt(3/2)+1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix2=sqrt(3/2)-1j*sqrt(dt*A)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));

tic;
Upre=U;
for t=1:nt
    
    if t==1
        
        F=U-dt*1/eps^2*laplacian9wM(laplacian9(U,N,h)-(2*U-6*U.^2+4*U.^3),N,h,mob(U))+A*dt*laplacian9(laplacian9(U,N,h),N,h);
        Fbar=fft2(F);
        Ubar=(Fbar./(invmatrix01))./(invmatrix02);
        U=real(ifft2(Ubar));
    else
        
        F=2*U-1/2*Upre+...
            2*(-dt*1/eps^2*laplacian9wM(laplacian9(U,N,h)-(2*U-6*U.^2+4*U.^3),N,h,mob(U))+A*dt*laplacian9(laplacian9(U,N,h),N,h))-...
        (-dt*1/eps^2*laplacian9wM(laplacian9(Upre,N,h)-(2*Upre-6*Upre.^2+4*Upre.^3),N,h,mob(Upre))+A*dt*laplacian9(laplacian9(Upre,N,h),N,h));
        Upre=U;
        Fbar=fft2(F);
        Ubar=(Fbar./(invmatrix1))./(invmatrix2);
        U=real(ifft2(Ubar));
        
    end
    if mod(t,2^15)==0
        save(['results/siechwm' num2str(t) 't' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps')
    end
end
toc;
save(['../results/siechwm' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps')