function [A,B,C,D,E,F,G,H,I]=getconsistvars(gamma,theta,theend)
N=size(gamma,1);
A=zeros(1,N);
B=zeros(1,N);
C=zeros(1,N);
D=zeros(1,N);
E=zeros(1,N);
F=zeros(1,N);
G=zeros(1,N);
H=zeros(1,N);
I=zeros(1,N);
A(1)=1/gamma(1,1);
B(1)=1/gamma(1,1)^2;
C(1)=0;
D(1)=1/2*1/gamma(1,1)^3;
E(1)=0;
F(1)=1/gamma(1,1)^3;
G(1)=0;
H(1)=0;
I(1)=0;
for i=2:N
    S=sum(gamma(i,1:N));
    theA=sum(A(1:i-1).*gamma(i,2:i));
    theB=sum(B(1:i-1).*gamma(i,2:i));
    theC=sum(C(1:i-1).*gamma(i,2:i));
    theD=sum(D(1:i-1).*gamma(i,2:i));
    theE=sum(E(1:i-1).*gamma(i,2:i));
    theF=sum(F(1:i-1).*gamma(i,2:i));
    theG=sum(G(1:i-1).*gamma(i,2:i));
    theH=sum(H(1:i-1).*gamma(i,2:i));
    theI=sum(I(1:i-1).*gamma(i,2:i));
    A(i)=(1+theA)/S;
    B(i)=(A(i)+theB)/S;
    C(i)=(sum(theta(i,2:i).*A(1:i-1))+theC)/S;
    %C(i)=(A(i-1)+theC)/S;
    D(i)=(A(i)^2/2+theD)/S;
    E(i)=(sum(theta(i,2:i).*A(1:i-1).^2)/2+theE)/S;
    F(i)=(B(i)+theF)/S;
    G(i)=(C(i)+theG)/S;
    H(i)=(sum(theta(i,2:i).*B(1:i-1))+theH)/S;
    I(i)=(sum(theta(i,2:i).*C(1:i-1))+theI)/S;
    
end

if theend==1
    A=A(N);
    B=B(N);
    C=C(N);
    D=D(N);
    E=E(N);
    F=F(N);
    G=G(N);
    H=H(N);
    I=I(N);
end

end