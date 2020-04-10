nt=2^23;
%N=1024;
% %nt=2^19;
N=ceil(sqrt(2)*1024);
% %nt=2^20;
%N=2048;
T=20;
eps=1;
[U,~,~,h]=initializebigoval(N,eps);
Uinit=U;
dt=T/nt;
[l,k]=meshgrid(0:(N-1));
invmatrix01=1+1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix02=1-1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix1=sqrt(3/2)+1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
invmatrix2=sqrt(3/2)-1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));

tic;
Upre=U;
for t=1:nt
    if t==1
        F=U+dt*1/eps^2*laplacian9(2*U-6*U.^2+4*U.^3,N,h);
        Fbar=fft2(F);
        Ubar=(Fbar./(invmatrix01))./(invmatrix02);
        U=real(ifft2(Ubar));
    else
        F=2*U-1/2*Upre+2*dt*1/eps^2*laplacian9(2*U-6*U.^2+4*U.^3,N,h)-dt*1/eps^2*laplacian9(2*Upre-6*Upre.^2+4*Upre.^3,N,h);
        Upre=U;
        Fbar=fft2(F);
        Ubar=(Fbar./(invmatrix1))./(invmatrix2);
        U=real(ifft2(Ubar));
    end
    if mod(t,2^22)==0
        save(['results/siech' num2str(nt) 'N' num2str(N) 't' num2str(t)],'T','dt','N','U','Upre','eps','t');
    end
end
toc;
save(['results/siech' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
