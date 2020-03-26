function y=applyoperator(x,MU,N,dt,eps,h,S)
    U=reshape(x, [N,N]);
    L=laplacian9(U,N,h);
    LW=laplacian9wM(eps*L+U,N,h,MU);
    y=S*U+dt*LW;
    y=y(:);
end