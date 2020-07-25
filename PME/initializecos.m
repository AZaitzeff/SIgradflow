function [u,x,h] = initializecos(N)
h=2/(N);
x=-1:h:(1-h);
u=cos(x*pi)+3/2;
end
