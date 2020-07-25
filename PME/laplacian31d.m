function L=laplacian31d(U,N,h)


    U=padarray(U,[0,2],'symmetric');
    L=(-2*U(3:N+2)+U(2:N+1)+U(4:N+3))/h^2;