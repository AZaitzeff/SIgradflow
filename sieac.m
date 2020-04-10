
%nt=2^22;
%N=4096;
nt=2^20;
N=2048;
T=20;
eps=1;
[U,X,Y,h]=initializebigcircle(N,eps);
Uinit=U;
dt=T/nt;
[l,k]=meshgrid(0:(N-1));
invmatrix0=1-dt/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix=3/2-dt/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
%invmatrix=1-dt/(2*h^2)*(-4976/1152+3/2*(cos(2*pi*l/N)+cos(2*pi*k/N))...
%    +3/40*(cos(4*pi*l/N)+cos(4*pi*k/N))+1/45*(cos(6*pi*l/N)+cos(6*pi*k/N))...
%    +2*(cos(2*pi*l/N)*cos(2*pi*k/N))+1/8*(cos(4*pi*l/N)*cos(4*pi*k/N))...
%    -1/2*(cos(4*pi*l/N)*cos(2*pi*k/N)+cos(2*pi*l/N)*cos(4*pi*k/N)));
tic;
Upre=U;
for t=1:nt
    if t==1
        F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
        Fbar=fft2(F);
        Ubar=Fbar./(invmatrix0);
        U=real(ifft2(Ubar));
    else
        
        F=2*U-1/2*Upre-2*dt*1/eps^2*(2*U-6*U.^2+4*U.^3)+dt*1/eps^2*(2*Upre-6*Upre.^2+4*Upre.^3);
        Upre=U;
        Fbar=fft2(F);
        Ubar=Fbar./(invmatrix);
        U=real(ifft2(Ubar));
    end
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    
    
end
toc;
save(['results/si2ac' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps')
