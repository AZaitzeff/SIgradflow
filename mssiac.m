
load('si3rd.mat')
fac=0;
T=20;
order=3;
eps=1;

nts=2.^(7:9);

 

for numt=1:3

    %N=ceil(512*2^((numt-2)/2));
    N=2048;
    [l,k]=meshgrid(0:(N-1));
    Uall=zeros(m+1,N,N);
    [Uinit,X,Y,h]=initializebigcircle(N,eps);
    nt=nts(numt);
    U=Uinit;
    dt=T/nt;
    tic;
for t=1:nt
    curU=U;
    Uall(1,:,:)=U;
    for step=1:m
        Utemp=squeeze(Uall(1,:,:));
        Un=Utemp*gamma(step,1)-dt*theta(step,1)*(2*(1-fac)*Utemp-6*Utemp.^2+4*Utemp.^3);

        for zm=2:step
            Utemp=squeeze(Uall(zm,:,:));
            Un=Un+Utemp*gamma(step,zm)-dt*theta(step,zm)*(2*(1-fac)*Utemp-6*Utemp.^2+4*Utemp.^3);
        end

        S=sum(gamma(step,:));
        invmatrix=S+2*dt*fac-dt/(h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
        Fbar=fft2(Un);
        Ubar=Fbar./(invmatrix);
        curU=real(ifft2(Ubar));
        
        Uall(step+1,:,:)=curU;

    end

    U=curU;

    

end

toc;

save(['results/multistepac' num2str(order) 's' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');

end