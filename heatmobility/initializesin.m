function [u,x,h] = initializesin(N)
h=1/(N);
x=(0+h/2):h:(1-h/2);
u=cos(x*(pi))+2;
end

