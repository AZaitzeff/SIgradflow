function L=laplacian5wM(U,N,h,M)
    L=dx(dx(U,N,h).*M,N,h)+dy(dy(U,N,h).*M,N,h);
    
end


function L=dx(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=1/2*(-U(3:N+2,2:N+1)+U(3:N+2,4:N+3))/h;
    
end

function L=dy(U,N,h)

    U=padarray(U,[2,2],'circular');
    L=1/2*(-U(2:N+1,3:N+2)+U(4:N+3,3:N+2))/h;
    
end