function UF=applyoperator(U,Us,dt,N,h,S)
    UF=S*U-dt*firstder51d(firstder51d(2*U,N,h).*Us,N,h);
end