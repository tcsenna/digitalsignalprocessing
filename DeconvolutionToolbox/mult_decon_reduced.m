function [Xk,red,J,J1,J2,perc]=mult_decon_reduced(x,Niter,Nline,epi,lambda,Nw,pcorr0)
%% 01/11/2016
[Ns, Nx] = size(x);
L        = Ns-Nw+1;
Xk       = x(1:L,:);
Xk       = Xk/norm(Xk,'fro');
J        = zeros(1,Niter);
J1       = zeros(1,Niter);
J2       = zeros(1,Niter);
Df       = zeros(L+Ns-1,Nx);
for ii=1:Nx
    Df(:,ii)=fft(x(:,ii),L+Ns-1);
end
aux_stop  = 0;
pcorr1    = pcorr0;
while aux_stop == 0
    idx3 = select_traces(x,pcorr1);
    if (length(unique(abs(idx3(:,2:5))))) == Nx
        aux_stop = 1;
    end
    pcorr1 = pcorr1 + 0.001;
end
red = size(idx3,1)/(Nx*(Nx-1)*2);
perc = size(idx3,1)/(Nx*(Nx-1)*2)*100;
[pcorr1-0.001 size(idx3,1) perc]
for kk=1:Niter
    
    %     kk/Niter
    
    Gq  = multiplyR3(Xk,Df,idx3,Nx,L,Ns); %%% A'Ax
    G   = Gq+lambda*Xk./((Xk.^2+epi^2).^(.5));
    Aux = G-sum(sum(Xk.*G))*Xk;
    Hk  = Aux/norm(Aux,'fro');
    tha = -1*pi/180;
    thb = 0;
    
    Xteste = cos(thb)*Xk+sin(thb)*Hk;
    fb     = 0.5*sum(sum(Xteste.*multiplyR3(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));
    phi    = 2-(1+5^.5)/2; % golden section search
    
    for jj=1:Nline
        
        thc=tha+phi*(thb-tha);
        thd=thb+phi*(tha-thb);
        
        Xteste=cos(thc)*Xk+sin(thc)*Hk;
        fc=0.5*sum(sum(Xteste.*multiplyR3(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));
        
        Xteste=cos(thd)*Xk+sin(thd)*Hk;
        fd=0.5*sum(sum(Xteste.*multiplyR3(Xteste,Df,idx3,Nx,L,Ns)))+lambda*sum(sum(((Xteste.^2+epi^2).^(.5))-epi));
        
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
    J1(kk) = 0.5*sum(sum(Xk.*multiplyR3(Xk,Df,idx3,Nx,L,Ns)));
    J2(kk) = +lambda*sum(sum(((Xk.^2+epi^2).^(.5))-epi));
    
end

end


