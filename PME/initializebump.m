function [u,x,h] = initializebump(N)
h=6/(N-1);
x=-3:h:(3);
u=x*0+1/100;
u(abs(x)<1)=exp(-1./(1-x(abs(x)<1).^2))+1/100;
end
