%matriz que cria full matrix A
%inputs são: dado completo, o tamanho da refletividade, e Nx sendo u número
%de traços que desejo usar
function[A]=criaA(x_awgn,L,Nx)
[Ns,Ntr]=size(x_awgn);
idx0     = zeros(Nx*(Nx-1)/2,Nx);
kk       = 1;
for ii=1:Nx
    for jj=ii+1:Nx
        idx0(kk,ii)   = jj;
        idx0(kk,jj)   = -ii;
        kk           = kk+1;
    end
end
[li,ci]=size(idx0);
%cada submatrizinha eh D

m=Ns+L-1;
p=L;
D=zeros(m,p);
A=zeros(m*li,p*ci);
%size(A)
for ll=1:li
    for cc=1:ci
        if (idx0(ll,cc)~=0)
            tr=abs(idx0(ll,cc)); %número do traço invocado
            s=sign(idx0(ll,cc)); %sinal do traço na matriz
            D=convmtx(x_awgn(:,tr),L);
            A(m*(ll-1)+1:m*ll,p*(cc-1)+1:p*cc)=s.*D;
        else
            A(m*(ll-1)+1:m*ll,p*(cc-1)+1:p*cc)=zeros(m,p);
        end
    end
end
