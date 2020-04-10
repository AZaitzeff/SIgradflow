%load('si2nd.mat')
load('si3rd.mat')
m=size(gamma,1);
fac=5;

n=2^13+1;
h=20/(n-1);
x=-10:h:10;
eps=1/4;
c=2;
T=5;
uinit=tanh((x+5)/eps);
trueu=tanh(((x+5)-c*T)/eps);
nts=[floor(702/m)];%,2^14,2^15];
%nts=[2^12,2^13,2^14,2^15,2^16];
N=numel(nts);
error1=zeros(2,N);

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
            utemp=U(1,:)*gamma(step,1)-dt*theta(step,1).*(c/eps-2*U(1,:)/eps^2-c/eps*U(1,:).^2+2/eps^2*U(1,:).^3-2*fac*U(1,:));
            for zm=2:step
                utemp=utemp+U(zm,:)*gamma(step,zm)...
                    -dt*theta(step,zm).*(c/eps-2*U(zm,:)/eps^2-c/eps*U(zm,:).^2+2/eps^2*U(zm,:).^3-2*fac*U(zm,:));
            end
            S=sum(gamma(step,:));
            utemp(3)=utemp(3)-4/3*dt/h^2+1/12*dt/h^2;
            utemp(4)=utemp(4)+1/12*dt/h^2;
            utemp(n-2)=utemp(n-2)+4/3*dt/h^2-1/12*dt/h^2;
            utemp(n-3)=utemp(n-3)-1/12*dt/h^2;

            B=zeros(n-4,5);
            B(:,3)=S+5/2*dt/h^2+2*fac*dt;
            B(:,2)=-4/3*dt/h^2;
            B(:,4)=-4/3*dt/h^2;
            B(:,1)=1/12*dt/h^2;
            B(:,5)=1/12*dt/h^2;
            F=utemp(3:n-2);
            A = spdiags(B,[-2,-1,0,1,2],n-4,n-4);
            ucur(3:n-2) = mldivide(A,F')';

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