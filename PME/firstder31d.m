function L=firstder31d(U,N,h)


    U=padarray(U,[2,0],'circular');
    L=(1/2*(-U(2:N+1)+U(4:N+3)))/h;