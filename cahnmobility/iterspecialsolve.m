function [U]=iterspecialsolve(U0,MU,dt,eps,N,h,S,f,tol)
[l,k]=meshgrid(0:(N-1));
%B=sqrt(S)-sqrt(dt*eps)/h^2*(-4+2*cos(2*pi*l/N)+2*cos(2*pi*k/N));
B=sqrt(S)-sqrt(dt*eps)/h^2*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));

maxiter=1000;
%A=sqrt(S)*speye(N^2)-sqrt(dt*eps)*laplacemtrxfast2orwM(N,h,MU(:));
A=sqrt(S)*speye(N^2)-sqrt(dt*eps)*laplacemtrxwM4ord(N,h,MU(:));
L = ichol(A,struct('michol','on'));
UB=U0-sqrt(eps*dt)*laplacian9(U0,N,h);
U=U0;
func=@(x) applyoperator(x,MU,N,dt,eps,h,S);
totiter=0;
tol1=max(abs(func(U(:))-f(:)))*.75;
for i=1:maxiter
    RHS=f-sqrt(S*dt*eps)*laplacian9wM(U,N,h,MU)-sqrt(S*dt*eps)*laplacian9(U,N,h);
    %RHS=f-sqrt(S*dt*eps)*laplacian5wM(U,N,h,MU)-sqrt(S*dt*eps)*laplacian5(U,N,h);
    %[x,flag,~,iter]=pcg(A,RHS(:),tol1/10,100,L,L',U(:));
    [x,iter]=zpcg(UB(:),A,RHS(:),tol1,100,L);
    totiter=totiter+iter;
    UB=reshape(x, [N,N]);
    UM1=real(ifft2(fft2(UB)./B));
    U=UM1;
    val=max(abs(func(U(:))-f(:)));
    if val<tol
        break
    end
    %val
    tol1=val*.75;
end