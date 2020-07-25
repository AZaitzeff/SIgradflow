function L=laplacian71d(U,N,h)


    U=padarray(U,[0,3],'symmetric');
    L=(-49/18*U(4:N+3)+3/2*(U(3:N+2)+U(5:N+4))-...
        3/20*(U(2:N+1)+U(6:N+5))+1/90*(U(1:N)+U(7:N+6)))/h^2;