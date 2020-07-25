
if order==2
    load('si2nd.mat')
else order==3
    load('si3rd.mat')
end
T=20;
fac=0;
eps=1;
order=3;
nts=2.^(7:11);


for numt=1:5
    N=2048;
    [l,k]=meshgrid(0:(N-1));
    Uall=zeros(m+1,N,N);
    [Uinit,X,Y,h]=initializebigoval(N,eps);
    nt=nts(numt);
    U=Uinit;
    dt=T/nt;
    tic;
for t=1:nt
    curU=U;
    Uall(1,:,:)=U;
    for step=1:m
        Utemp=squeeze(Uall(1,:,:));
        Un=Utemp*gamma(step,1)+dt*theta(step,1)*laplacian9(2*(1-fac)*Utemp-6*Utemp.^2+4*Utemp.^3,N,h);

        S=sum(gamma(step,:));
        
        for zm=2:step
            Utemp=squeeze(Uall(zm,:,:));
            Un=Un+Utemp*gamma(step,zm)+dt*theta(step,zm)*laplacian9(2*(1-fac)*Utemp-6*Utemp.^2+4*Utemp.^3,N,h);
        end
        
        
        invmatrix1=sqrt(S)+1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
        invmatrix2=sqrt(S)-1j*sqrt(dt)/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
        
        Fbar=fft2(Un);
        Ubar=(Fbar./(invmatrix1))./(invmatrix2);
        curU=real(ifft2(Ubar));
        %val=normdiff(curU,Un,N,dt,eps,h,S,force)
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
end
toc;
save(['results/multistepch' num2str(order) 's' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
%sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum((Uinit(:)-trueU(:)).^2))
%error1(1,numt)=sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum(trueU(:).^2));
end