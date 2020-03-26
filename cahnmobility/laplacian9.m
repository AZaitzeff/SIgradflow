function L=laplacian9(U,N,h)


    U=padarray(U,[2,2],'circular');
    L=(-5*U(3:N+2,3:N+2)+4/3*(U(2:N+1,3:N+2)+U(4:N+3,3:N+2)+U(3:N+2,2:N+1)+U(3:N+2,4:N+3))-...
        1/12*(U(1:N,3:N+2)+U(5:N+4,3:N+2)+U(3:N+2,1:N)+U(3:N+2,5:N+4)))/h^2;