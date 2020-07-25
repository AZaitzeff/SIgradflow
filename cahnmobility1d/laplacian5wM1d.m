function L=laplacian5wM1d(U,N,h,M)
    L=dx(dx(U,N,h).*M,N,h);
    
end


function L=dx(U,N,h)

    U=padarray(U,[0,2],'circular');
    L=(2/3*(-U(2:N+1)+U(4:N+3))-1/12*(-U(1:N)+U(5:N+4)))/h;
    
end