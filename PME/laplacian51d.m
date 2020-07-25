function L=laplacian51d(U,N,h)


    U=padarray(U,[0,2],'symmetric');
    L=(-5/2*U(3:N+2)+4/3*(U(2:N+1)+U(4:N+3))-...
        1/12*(U(1:N)+U(5:N+4)))/h^2;