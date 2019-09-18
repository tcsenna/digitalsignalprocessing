Nx=4;
B=[1 2;3 4;5 6];
[m,p]=size(B);
idx0     = zeros(Nx*(Nx-1)/2,2);
[li,ci]=size(idx0);
A=zeros(li*m,Nx*p);
kk       = 1;
for ii=1:Nx
    for jj=ii+1:Nx
        idx0(kk,1)   = jj;
        idx0(kk,2)   = -ii;
        kk           = kk+1;
    end
end
la=1;
ca=1;
for kk=1:li
    if (kk<=(li-1)) && (idx0(kk,2)== idx0(kk+1,2))%vou preencher até a linha kk-1 de idx0  
        A(m*(kk-1)+1:kk*m,p*(kk-1)+1:kk*p)=B
        A(m*(kk-1)+1:
    
    
    
    
    
    
    
    
    
    end
    
end
