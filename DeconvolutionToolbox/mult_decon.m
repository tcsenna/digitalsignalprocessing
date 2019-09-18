function [Xk,J,J1,J2]=mult_decon(x,Niter,Nline,epi,lambda,Nw)
%% 08/11/2016
[Ns, Nx] = size(x);
L        = Ns-Nw+1;
Xk       = x(1:L,:);
Xk       = Xk/norm(Xk,'fro');
% J        = zeros(1,Niter);
% J1       = zeros(1,Niter);
% J2       = zeros(1,Niter);
idx0     = zeros(Nx*(Nx-1)/2,2);
kk       = 1;
for ii=1:Nx
    for jj=ii+1:Nx
        idx0(kk,1)   = jj;
        idx0(kk,2)   = -ii;
        kk           = kk+1;
    end
end
%acima criei a matriz que representa todos os pares possíveis de traços
% figure(100), imagesc(abs(idx0))
Df = zeros(L+Ns-1,Nx);
for ii=1:Nx
    Df(:,ii)=fft(x(:,ii),L+Ns-1);
end
%acima a matriz Df tem o tamanho da convolucao entre a refletividade e o
%traço e os elementos sao a fft de cada traço desse tamanho (convolucao entre refletividade e traço) 
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
size(idx3,1)
for kk=1:Niter
    
%     kk/Niter
    
    Gq  = multiplyR(Xk,Df,idx3,Nx,L,Ns); %%% A'Ax
    G   = Gq+lambda*Xk./((Xk.^2+epi^2).^(.5));
    %
    Aux = G-sum(sum(Xk.*G))*Xk;     
    Hk  = Aux/norm(Aux,'fro');  
    tha = -1*pi/180;
    thb = 0;    
    
    Xteste = cos(thb)*Xk+sin(thb)*Hk;
    fb     = 0.5*sum(sum(Xteste.*multiplyR(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));
    phi    = 2-(1+5^.5)/2; % golden section search

    for jj=1:Nline
        
        thc=tha+phi*(thb-tha);
        thd=thb+phi*(tha-thb);
       
        Xteste=cos(thc)*Xk+sin(thc)*Hk;
        fc=0.5*sum(sum(Xteste.*multiplyR(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));        
        
        Xteste=cos(thd)*Xk+sin(thd)*Hk;        
        fd=0.5*sum(sum(Xteste.*multiplyR(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));        

        if (fc > fd)
            
            tha=thc;  
        else
            
            thb=thd;
            fb=fd;                    
            
        end
        
    end

    if(thb==0)
        break;
    end
    
    Xk     = cos(thb)*Xk+sin(thb)*Hk;
    J(kk)  = fb;
    J1(kk) = 0.5*sum(sum(Xk.*multiplyR(Xk,Df,idx3,Nx,L,Ns)));
    J2(kk) = +lambda*sum(sum(((Xk.^2+epi^2).^(.5))-epi));

end

end


