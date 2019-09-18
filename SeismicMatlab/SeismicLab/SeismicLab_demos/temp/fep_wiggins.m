function [e_wiggins,h_wiggins,custo]=fep_wiggins(entrada_temp,h,N_iter,mu,sinal,E)

%%% Filtro de erro de predi��o (FEP)
%%% Estrutura: Linear
%%% Crit�rio: Wiggins (Kurtosis)
%%% entrada=vetor de entrada do FEP (deve ser um vetor linha)
%%% N=n�mero de amostras do filtro de predi��o (coeficientes n�o nulos)
%%% N_iter=n�mero de itera��es (6 a 10 itera��es)
%%% E=% de branqueamento de R


mincorr=1;
maxcorr=size(entrada_temp,1);
entrada=entrada_temp(mincorr:maxcorr,:);

% h=zeros(1,N);
% h(floor(N/2))=1;
% h=h./sqrt(sum(h.^2));

N=size(h,2);

I=eye(N);

for ll=1:N_iter
    
    ll/N_iter
    
    RR=zeros(N,N);
    pp=zeros(N,1);
    vv(ll)=0;
    
    for jj=1:size(entrada,2)
        
        kk=1;
        for ii=1:N
            U(:,kk)=filter(I(ii,:),1,entrada(:,jj));
            kk=kk+1;
        end
        
        y=filter(h,1,entrada(:,jj));
        y3=y.^3;
        v(ll)=sum(y.^4)./((sum(y.^2)).^2);
        s=sum(y.^2);
        
        R=(v(ll)/s)*U'*U;
        R=R+R(1,1)*E/100;
        
        p=(1/(s.^2))*U'*y3;
        RR=RR+R;
        pp=pp+p;
        vv(ll)=vv(ll)+v(ll);
        
    end
    
    h_temp=pinv(RR)*pp;
    h_temp=h_temp./sqrt(sum(h_temp.^2));
    h=h+sinal*mu*h_temp';
    h=h./sqrt(sum(h.^2));
    
end

figure
plot(vv)

custo=vv(N_iter);

h_wiggins=h;
e_wiggins=filter(h_wiggins,1,entrada_temp);
end




