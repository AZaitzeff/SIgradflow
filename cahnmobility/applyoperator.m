function y=applyoperator(x,MU,N,dt,eps,h,S)
    U=reshape(x, [N,N]);
    L=laplacian5(U,N,h);
    LW=laplacian5wM(eps*L,N,h,MU);
    y=S*U+dt*LW;
    y=y(:);
end