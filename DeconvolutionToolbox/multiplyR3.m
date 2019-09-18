function [Xout]=multiplyR3(Xk,Df,idx3,Nx,L,Ns)
%% 01/11/2016

Xf   = zeros(L+Ns-1,Nx);
Xout = zeros(L,Nx);

for jj=1:Nx
    Xf(:,jj) = fft(Xk(:,jj),L+Ns-1);
end

for ll=1:size(idx3,1)
    kk         = idx3(ll,1);
    jj         = idx3(ll,2);
    ii         = idx3(ll,3);
    mm         = idx3(ll,4);
    nn         = idx3(ll,5);
    aux        = sign(ii)*sign(jj)*ifft(conj(Df(:,abs(ii))).*Df(:,abs(jj)).*Xf(:,mm));
    Xout(:,nn) = Xout(:,nn)+aux(1:L,1);
end

% for nn=1:Nx;
%     NN = sum(idx3(:,5)==nn)/2;
%     a  = ((Nx-1)/NN);
%     Xout(:,nn) = a*Xout(:,nn);
% end

end