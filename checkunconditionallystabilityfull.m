function [cond,flag]=checkunconditionallystabilityfull(gamma,theta,z)
M=size(gamma,1);
gammat=zeros(M,M);
S=zeros(M,M);


gammat(M,:)= gamma(M,:);
S(M,:)=cumsum(gammat(M,:));
for m=(M-1):-1:1
   for i=1:m
       val=0;
       for j=(m+1):M
           val=val+gammat(j,i)*S(j,m)/S(j,j);
       end
      gammat(m,i)=gamma(m,i)-val-z*theta(m,i);
   end
   S(m,1:m)=cumsum(gammat(m,1:m));
end

cond=diag(S);
if any(cond<0)
    flag=0;
else
    flag=1;
end