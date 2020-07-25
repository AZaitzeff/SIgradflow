function L=firstdermtrxfast(N,h)
m=N;
%L = sparse(m,m);
coefs=[0,1/12,-2/3,2/3,-1/12];
rows=zeros(1,5*m);
cols=zeros(1,5*m);
vals=zeros(1,5*m);

indices=1:m;
rows(1:m)=indices;
cols(1:m)=indices;
vals(1:m)=0;

for i=2:5
    rows((i-1)*m+1:i*m)=indices;
    if i==2
        newindices=indices-2;
        change=1;
        newindices(change)=newindices(change)+N;
        change=2;
        newindices(change)=newindices(change)+N;
    elseif i==3
        newindices=indices-1;
        change=1;
        newindices(change)=newindices(change)+N;
    elseif i==4
        newindices=indices+1;
        change=N;
        newindices(change)=newindices(change)-N;
    elseif i==5
        newindices=indices+2;
        change=N;
        newindices(change)=newindices(change)-N;
        change=(N-1);
        newindices(change)=newindices(change)-N;
    end
    
    cols((i-1)*m+1:i*m)=newindices;
    vals((i-1)*m+1:i*m)=coefs(i)/h;
    
end




L = sparse(rows,cols,vals,m,m);