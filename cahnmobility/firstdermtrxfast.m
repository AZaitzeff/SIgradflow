function [Lx,Ly]=firstdermtrxfast(N,h)
m=N^2;
%L = sparse(m,m);
coefs=[1/12,-2/3,2/3,-1/12,1/12,-2/3,2/3,-1/12];
rows=zeros(1,8*m);
cols=zeros(1,8*m);
vals=zeros(1,8*m);
indices=1:m;

for i=1:8
    rows((i-1)*m+1:i*m)=indices;
    if i==1
        newindices=indices-2;
        change=1:N:m;
        newindices(change)=newindices(change)+N;
        change=2:N:m;
        newindices(change)=newindices(change)+N;
    elseif i==2
        newindices=indices-1;
        change=1:N:m;
        newindices(change)=newindices(change)+N;
    elseif i==3
        newindices=indices+1;
        change=N:N:m;
        newindices(change)=newindices(change)-N;
    elseif i==4
        newindices=indices+2;
        change=N:N:m;
        newindices(change)=newindices(change)-N;
        change=(N-1):N:m;
        newindices(change)=newindices(change)-N;
    elseif i==5
        newindices=indices-2*N;
        change=1:N;
        newindices(change)=newindices(change)+N^2;
        change=(N+1):2*N;
        newindices(change)=newindices(change)+N^2;
        
    elseif i==6
        newindices=indices-N;
        change=1:N;
        newindices(change)=newindices(change)+N^2;
        
    elseif i==7
        newindices=indices+N;
        change=(N*(N-1)+1):m;
        newindices(change)=newindices(change)-N^2;
        
    elseif i==8
        newindices=indices+2*N;
        change=(N*(N-1)+1):m;
        newindices(change)=newindices(change)-N^2;
        change=(N*(N-2)+1):N*(N-1);
        newindices(change)=newindices(change)-N^2;
    end
    
    cols((i-1)*m+1:i*m)=newindices;
    vals((i-1)*m+1:i*m)=coefs(i)/h;
    
end

Ly = sparse(rows(1:4*m),cols(1:4*m),vals(1:4*m),m,m);
Lx = sparse(rows(4*m+1:8*m),cols(4*m+1:8*m),vals(4*m+1:8*m),m,m);