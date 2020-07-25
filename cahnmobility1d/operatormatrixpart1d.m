function L=operatormatrixpart1d(MU,N,dt,eps,h)
    DD=laplacemtrxfast(N,h);
    D=firstdermtrxfast(N,h);
    MUdaig=spdiags(MU',[0],N,N);
    %LW=laplacian9wM(eps*L,N,h,MU);
    L=dt*eps*D*MUdaig*D*DD;
end