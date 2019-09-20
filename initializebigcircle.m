function [U,X,Y,h]=initializebigcircle(N,eps)
h=20/(N);
x=-10:h:10-h;
[X,Y] = meshgrid(x);
U=1./(1+exp(-1/eps^2*(7.5-sqrt(X.^2+Y.^2))));