function [u,x,h] = initializesin(N)
h=1/N;
x=0:h:1-h;
u=sin(x*(2*pi))+2;
end

