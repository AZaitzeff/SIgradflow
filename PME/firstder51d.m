function L=firstder51d(U,N,h)


    U=padarray(U,[2,0],'circular');
    L=(2/3*(-U(2:N+1)+U(4:N+3))+...
        1/12*(U(1:N)-U(5:N+4)))/h;