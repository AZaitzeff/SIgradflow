function [U,x,h]=initialize1d(N,eps)
h=1/(N);
x=-.5:h:.5-h;
U=tanh(1/(10*sqrt(eps))*cos(2*pi*x));