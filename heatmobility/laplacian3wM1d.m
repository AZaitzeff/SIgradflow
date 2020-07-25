function L=laplacian3wM1d(U,N,h,M)
    Ms=padarray(M,[0,1],'symmetric');
    Ms12=(Ms(1:(N+1))+Ms(2:(N+2)))/2;
    Us=padarray(U,[0,1],'symmetric');
    L=(Us(1:N).*Ms12(1:N)-Us(2:(N+1)).*(Ms12(1:N)+Ms12(2:(N+1)))+Us(3:(N+2)).*Ms12(2:(N+1)))/h^2;
    
    
end

