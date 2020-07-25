function [u,x,h] = initializegauss(N,sig)
h=6/(N-1);
x=-3:h:(3);
u=exp(-1/2*(x/sig).^2)/(sig*sqrt(2*pi));
end
