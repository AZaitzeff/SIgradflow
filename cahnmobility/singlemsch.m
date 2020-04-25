function [UF,RHS]=singlemsch(U0,gamma,theta,m,MU,dt,eps,N,h)
%A2=operatormatrixpart(MU,N,dt,eps,h);
tol1=1e-4;
E2=@(U) dt*laplacian5wM((U.^2-1).*U,N,h,MU);
Uall=zeros(m,N,N);
UF=U0;
for step=1:m
    RHS=U0*gamma(step,1)+E2(U0)*theta(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:,:));
        RHS=RHS+curU*gamma(step,zm)+E2(curU)*theta(step,zm);
    end
    S=sum(gamma(step,:));
    [UF]=iterspecialsolve(UF,MU,dt,eps,N,h,S,RHS,tol1);
    %RHSs=1/S*dt*eps*laplacian5wM(laplacian5(RHS,N,h),N,h,MU);
    %func=@(x) applyoperatorbreak(x,MU,N,dt,eps,h,S);
    
    
    %Us=RHS-UF;
    %Us=laplacian5(Us,N,h);
    %tic
    %[x]=zcg(Us(:),func,RHSs(:),tol1,3000);
    %toc
    %max(abs(func(x)-RHSs(:)))
    %[x,flag,res,~,fill]=pcg(func,RHSs(:),tol1,1000,[],[]);
    %flag
    %max(abs(func(x)-RHSs(:)))
    %plot(fill)
    %flag
    
    %A=S*speye(N^2)+A2;
    %L1 = ichol(A, struct('type','ict','diagcomp',65));
    %tic
    %[L,Up] = ilu(A);
    %[x,flag,~,iter] = gmres(A,RHS(:),40,tol1,3000,L,Up,UF(:));
    %toc
    %iter
    %x=A\RHS(:);
    %if flag>0
    %   flag
    %   N
    %   dt
    %   'warning' 
    %end
    %UF=reshape(x, [N,N]);
    func=@(x) applyoperator(x,MU,N,dt,eps,h,S);
    max(abs(func(UF(:))-RHS(:)))
    %max(abs(func(UF1(:))-RHS(:)))
    %max(abs(func(x)-RHS(:)))
    %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    Uall(step,:,:)=UF;
end


