function L=laplacemtrxfast2orwM(N,h,MU)
m=N^2;
%L = sparse(m,m);
rows=zeros(5*m,1);
cols=zeros(5*m,1);
vals=zeros(5*m,1);

indices=1:m;
rows(1:m)=indices;
cols(1:m)=indices;

for i=2:5
    rows((i-1)*m+1:i*m)=indices;
    if i==2
        newindices=indices-2;
        change=1:N:m;
        newindices(change)=newindices(change)+N;
        change=2:N:m;
        newindices(change)=newindices(change)+N;

        newindicesM=indices-1;
        change=1:N:m;
        newindicesM(change)=newindicesM(change)+N;
    elseif i==3
        newindicesM=indices+1;
        change=N:N:m;
        newindicesM(change)=newindicesM(change)-N;

        newindices=indices+2;
        change=N:N:m;
        newindices(change)=newindices(change)-N;
        change=(N-1):N:m;
        newindices(change)=newindices(change)-N;
    elseif i==4
        newindices=indices-2*N;
        change=1:N;
        newindices(change)=newindices(change)+N^2;
        change=(N+1):2*N;
        newindices(change)=newindices(change)+N^2;
        

        newindicesM=indices-N;
        change=1:N;
        newindicesM(change)=newindicesM(change)+N^2;
        
    elseif i==5
        newindicesM=indices+N;
        change=(N*(N-1)+1):m;
        newindicesM(change)=newindicesM(change)-N^2;
        

        newindices=indices+2*N;
        change=(N*(N-1)+1):m;
        newindices(change)=newindices(change)-N^2;
        change=(N*(N-2)+1):N*(N-1);
        newindices(change)=newindices(change)-N^2;
    end
    
    cols((i-1)*m+1:i*m)=newindices;
    vals((i-1)*m+1:i*m)=MU(newindicesM)/(4*h^2);
  
    vals(1:m)=vals(1:m)-MU(newindicesM)/(4*h^2);
    
end




L = sparse(rows,cols,vals,m,m);