function [UF]=singlemsch(U0,gamma,m,MU,dt,eps,N,h,C)

tol1=1e-8;

Uall=zeros(m,N,N);
UF=U0;
for step=1:m
    RHS=U0*gamma(step,1);
    for zm=2:step
        curU=squeeze(Uall(zm-1,:,:));
        RHS=RHS+curU*gamma(step,zm);
    end
    S=sum(gamma(step,:));
    if C>0
        [UF]=linearsolvefft(UF,MU,RHS,S,dt,eps,N,h,C);
        %func=@(x) applyoperator(x,MU,N,dt,eps,h,S);
        %max(abs(func(UF(:))-RHS(:)))
    else
        func=@(x) applyoperator(x,MU,N,dt,eps,h,S);
        [x,flag,~,~]=pcg(func,RHS(:),tol1,3000,[],[],UF(:));
        %max(abs(func(x)-RHS(:)))
        %iter
        if flag>0
           flag
           N
           dt
           'warning' 
        end
        UF=reshape(x, [N,N]);
    end
    %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    UF(UF>1)=1;
    UF(UF<-1)=-1;
    Uall(step,:,:)=UF;
end


