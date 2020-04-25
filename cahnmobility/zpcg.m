function [x,i]=zpcg(x,A,b,tol,maxiter,L)
r=b-A*x;
rb=L'\(L\(r));
p=rb;
for i=1:maxiter
    Ap=A*p;
    alpha=-rb'*r/(p'*Ap);
    x=x-alpha*p;
    rk1=r+alpha*Ap;
    %max(abs(rk1))<tol
    rk1b=L'\(L\(rk1));
    if max(abs(rk1))<tol
        break
    end
    p=rk1b+(rk1b'*rk1)/(rb'*r)*p;
    r=rk1;
    rb=rk1b;
    
    
end
%i