function [U,X,Y,h]=initializerandom(N)
rng(101)
h=1/(N);
x=0:h:1-h;
[X,Y] = meshgrid(x);
U=rand(N,N)*.1-.45;