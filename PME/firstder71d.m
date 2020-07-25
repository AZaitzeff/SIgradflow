function L=firstder71d(U,N,h)


    U=padarray(U,[3,0],'circular');
    L=(3/4*(-U(3:N+2)+U(5:N+4))+...
        3/20*(U(2:N+1)-U(6:N+5))+1/60*(-U(1:N)+U(7:N+6)))/h;