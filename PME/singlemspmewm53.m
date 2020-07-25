function [UF]=singlemspmewm53(U0,gamma,m,MU,dt,N,h)


tol=1e-8;
maxiter=1000;

Uall=zeros(m,N);
UF=U0;
for step=1:m
    RHS=U0*gamma(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:));
        RHS=RHS+curU*gamma(step,zm);
    end
    S=sum(gamma(step,:));
    

    
    for i=1:maxiter
        val=RHS-(S*UF-5/2*dt*laplacian3wM1d(nthroot(UF,3).^2,N,h,MU));
        
        if sqrt(sum(abs(val).^2))<tol
            break
        end
        [delUF,~]=onedsolvematrix(MU,1./(nthroot(UF,3)),h,S,5/3*dt,val,N);
        %[delUF,~]=onedfbsolve53matrix(MU,1./nthroot(UF,3),h,S,5/3*dt,val,N);
        UF=UF+delUF;
        
    end
    Uall(step,:)=UF;
end


