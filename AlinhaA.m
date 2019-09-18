%este programa calcula a matriz A, autovalores e autovetores de A'A, os
%dois menores autovalores de A'A e suas posições na matriz
%Ntr é o número de traços do dado que vc quer usar
%e_val1 e e_val2 sã0 o primeiro menor e segundo menor eigenvalues
[e_val,e_vec,e_val1,e_val2]=AlinhaA(x_awgn, Ns, Ntr, L)

A=zeros(((Ntr(Ntr-1))/2)*(Ns+L-1),Ntr*L);
kk       = 1;
for ii=1:Ntr
    for jj=ii+1:Ntr
        idx0(kk,1)   = jj;
        idx0(kk,2)   = -ii;
        kk           = kk+1;
    end
end
%linha de A e coluna de A
la=1;
for 
    
ll = 1;
for kk=1:(Nx*(Nx-1)/2)
    if ( sign(idx0(kk,1))~=0 ) && ( sign(idx0(kk,2))~=0 )
        idx3(ll,:)   = [kk  idx0(kk,1)  idx0(kk,1) abs(idx0(kk,2)) abs(idx0(kk,2))];
        idx3(ll+1,:) = [kk  idx0(kk,2)  idx0(kk,2) abs(idx0(kk,1)) abs(idx0(kk,1))];
        idx3(ll+2,:) = [kk  idx0(kk,1)  idx0(kk,2) abs(idx0(kk,2)) abs(idx0(kk,1))];
        idx3(ll+3,:) = [kk  idx0(kk,2)  idx0(kk,1) abs(idx0(kk,1)) abs(idx0(kk,2))];
        ll           = ll+4;
    end
end