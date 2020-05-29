function L=laplacemtrxwM4ord(N,h,MU)
%L = sparse(m,m);
[Lx,Ly]=firstdermtrxfast(N,h);
%[Lx,Ly]=firstdermtrxfast2or(N,h);
MUdaig=spdiags(MU,0,N^2,N^2);
L=(Lx*MUdaig*Lx+Ly*MUdaig*Ly);