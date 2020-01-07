function [uxf,uyf]=singlemultistep(ux0,uy0,gamma,theta,m,MU,k)

ux=zeros(1,m+1);
uy=zeros(1,m+1);
ux(1)=ux0;
uy(1)=uy0;
for j=2:m+1
    ux(j)=(sum(gamma(j-1,1:j-1).*ux(1:j-1))-0*MU*k*sum(ux(1:j-1).*theta(j-1,1:j-1)))/(sum(gamma(j-1,1:j-1))+MU*k);
    uy(j)=(sum(gamma(j-1,1:j-1).*uy(1:j-1))+MU*k*sum(uy(1:j-1).*theta(j-1,1:j-1)))/(sum(gamma(j-1,1:j-1))-0*MU*k);
end
uxf=ux(m+1);
uyf=uy(m+1);