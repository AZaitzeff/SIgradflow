function [U]=onednonconstantlinear(Us,h,S,k,f,n)
U=zeros(1,n);
c=zeros(1,n-1);
d=zeros(1,n);
delta=k/h^2;
c(1)=-delta*(Us(1)+Us(2))/2/(S+delta*(Us(1)+Us(2))/2);
d(1)=f(1)/(S+delta*(Us(1)+Us(2))/2);
for i=2:(n-1)
   c(i)= -delta*(Us(i)+Us(i+1))/2/((S+delta*(Us(i)+Us(i+1)/2+Us(i-1)/2))+delta*(Us(i-1)+Us(i))/2*c(i-1));
   d(i)= (f(i)+delta*(Us(i-1)+Us(i))/2*d(i-1))/((S+delta*(Us(i)+Us(i+1)/2+Us(i-1)/2))+delta*(Us(i-1)+Us(i))/2*c(i-1));
end

d(n)= (f(n)+delta*(Us(n-1)+Us(n))/2*d(n-1))/((S+delta*(Us(n)+Us(n-1))/2)+delta*(Us(n-1)+Us(n))/2*c(n-1));

U(n)=d(n);

for i=(n-1):-1:1
    U(i)=d(i)-c(i)*U(i+1);
    
end