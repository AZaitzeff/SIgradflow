function [U,X,Y,h]=initializebigoval(N,eps,theta)
h=1/(N);
x=0:h:1-h;
[X,Y] = meshgrid(x);
U=tanh(1/(10*sqrt(eps))*(1-(9*((X-.5)*cos(theta)+(Y-.5)*sin(theta)).^2+25*(-(X-.5)*sin(theta)+(Y-.5)*cos(theta)).^2)));