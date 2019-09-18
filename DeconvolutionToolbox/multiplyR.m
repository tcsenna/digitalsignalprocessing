function [Xout]=multiplyR(Xk,Df,idx3,Nx,L,Ns)
%% 08/11/2016

Xf   = zeros(L+Ns-1,Nx);
Xout = zeros(L,Nx);
%xk é uma matriz com as L primeiras linhas da matriz de dados x_awgn, normalizada
%pela norma de frobenius
%aqui eu calculo a fft de cada uma dessas colunas (traços) com L+Ns-1
%elementos, que é o tamanho da convolucao entre a refletividade e cada
%traço
for jj=1:Nx
    Xf(:,jj) = fft(Xk(:,jj),L+Ns-1);
end

for ll=1:size(idx3,1)
    jj         = idx3(ll,2);
    ii         = idx3(ll,3);
    mm         = idx3(ll,4);
    nn         = idx3(ll,5);
    aux        = sign(ii)*sign(jj)*ifft(conj(Df(:,abs(ii))).*Df(:,abs(jj)).*Xf(:,mm));
    Xout(:,nn) = Xout(:,nn)+aux(1:L,1);
end

end