function [U,X,Y,h]=initializebigoval(N,eps,theta)
h=1/(N);
x=-.5:h:.5-h;
[X,Y] = meshgrid(x);
U=tanh(1/(10*sqrt(eps))*(1-(9*((X)*cos(theta)+(Y)*sin(theta)).^2+25*(-(X)*sin(theta)+(Y)*cos(theta)).^2)));