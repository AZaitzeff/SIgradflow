function [u,x,h] = initializebarenblatt53(N,t)
h=8/(N-1);
x=-4:h:(4);
gam=5/3;
C=(sqrt((gam-1)/(2*pi*gam*(gam+1)))*gamma(3/2+1/(gam-1))/gamma(gam/(gam-1)))^((gam-1)/(gam+1));
u=(C^2-3/40*x.^2/(t)^(3/4));
u(u<0)=0;
u=1/(t)^(3/8)*u.^(3/2);
end