function L=laplacian3wMfb1d(U,N,h,M)

    Ms=padarray(M,[0,1],'symmetric');
    Us=padarray(U,[0,1],'symmetric');
    L=(Us(1:N).*Ms(2:(N+1))-Us(2:(N+1)).*(Ms(3:(N+2))+Ms(2:(N+1)))+Us(3:(N+2)).*Ms(3:(N+2)))/h^2;
    
    
end