function [U,A]=onedsolvematrix(Us,Un,h,S,k,f,n)

B=zeros(n,1);
B(2:n-1,:)=-k*(Us(2:n-1)+Us(1:n-2)/2+Us(3:n)/2)/h^2;
B(1)=-k*(Us(1)/2+Us(2)/2)/h^2;
B(n)=-k*(Us(n)/2+Us(n-1)/2)/h^2;

C=zeros(n,1);
C(2:n,:)=k*(Us(1:n-1)+Us(2:n))/(2*h^2);
A=zeros(n,1);
A(1:n-1,:)=k*(Us(1:n-1)+Us(2:n))/(2*h^2);

D = spdiags([A,B,C],[-1,0,1],n,n);

Undiag = spdiags(Un',0,n,n);
A=S*speye(n)-D*Undiag;
U=A\(f');
U=U';
