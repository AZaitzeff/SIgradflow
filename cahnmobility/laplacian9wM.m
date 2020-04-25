function L=laplacian9wM(U,N,h,M)
    %L=dx(U,N,h)+dy(U,N,h);
    %L=bdx(fdx(U,N,h).*M,N,h)+bdy(fdy(U,N,h).*M,N,h);
    L=dx(dx(U,N,h).*M,N,h)+dy(dy(U,N,h).*M,N,h);
    %L=dx(U,N,h).*dx(M,N,h)+dy(U,N,h).*dy(M,N,h)+laplacian9(U,N,h).*M;
    
end


function L=dx(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(2/3*(-U(3:N+2,2:N+1)+U(3:N+2,4:N+3))-...
        1/12*(-U(3:N+2,1:N)+U(3:N+2,5:N+4)))/h;
    
end

function L=dy(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(2/3*(-U(2:N+1,3:N+2)+U(4:N+3,3:N+2))-...
        1/12*(-U(1:N,3:N+2)+U(5:N+4,3:N+2)))/h;
    
end


function L=fdx(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(-3/2*U(3:N+2,3:N+2)+2*U(3:N+2,4:N+3)-1/2*U(3:N+2,5:N+4))/h;
    
end

function L=fdy(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(-3/2*U(3:N+2,3:N+2)+2*U(4:N+3,3:N+2)-1/2*U(5:N+4,3:N+2))/h;
    
end

function L=bdx(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(3/2*U(3:N+2,3:N+2)-2*U(3:N+2,2:N+1)+1/2*U(3:N+2,1:N))/h;
    
end

function L=bdy(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=(3/2*U(3:N+2,3:N+2)-2*U(2:N+1,3:N+2)+1/2*U(1:N,3:N+2))/h;
    
end