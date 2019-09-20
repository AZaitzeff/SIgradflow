addpath('../findgamma/')
[gamma,m]=getgamma(2,0);

n=2^13+1;
h=20/(n-1);
x=-10:h:10;
eps=1/4;
c=2;
T=5;
uinit=tanh((x+5)/eps);
trueu=tanh(((x+5)-c*T)/eps);
%nts=[2^7,2^8,2^9,2^10,2^11,2^12,2^13];%,2^14,2^15];
nts=[2^13,2^14,2^15,2^16];
N=numel(nts);
error1=zeros(2,N);
tol=1e-12;
U=zeros(m+1,n);
maxiter=100;
for numt=1:N
    nt=nts(numt);
    dt=T/nt;
    u=uinit;
    tic;
    for t=1:nt
        U(1,:)=u;
        ucur=u;
        for step=1:m
            utemp=U(1,:)*gamma(step,1);
            for zm=2:step
                utemp=utemp+U(zm,:)*gamma(step,zm);
            end
            S=sum(gamma(step,:));
            
            for z =1:maxiter
                ulast=ucur;
                B=zeros(n-4,5);
                B(:,3)=5/2*dt/h^2+S+dt*(-2/eps^2-2*c/eps*ucur(3:n-2)+6/eps^2*ucur(3:n-2).^2);
                B(:,2)=-4/3*dt/h^2;
                B(:,4)=-4/3*dt/h^2;
                B(:,1)=1/12*dt/h^2;
                B(:,5)=1/12*dt/h^2;
                F=-dt/h^2*(-5/2*ucur(3:n-2)+4/3*(ucur(2:n-3)+ucur(4:n-1))-1/12*(ucur(1:n-4)+ucur(5:n)))...
                    +S*ucur(3:n-2)+dt*(c/eps-2*ucur(3:n-2)/eps^2-c/eps*ucur(3:n-2).^2+2/eps^2*ucur(3:n-2).^3)-utemp(3:n-2);
                A = spdiags(B,[-2,-1,0,1,2],n-4,n-4);
                ucur(3:n-2) = ucur(3:n-2)- mldivide(A,F')';

                if sum(abs(ulast-ucur))<tol
                    break
                end
            end
            U(step+1,:)=ucur;
        end
    u=U(m+1,:);
    end
    toc;
    error1(1,numt)=sqrt(h*sum((trueu-u).^2));
end


for i=2:N
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end