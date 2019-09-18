function [e,melhor_h_g]=fep_l1l2_pso(entrada,N,N_iter,Np,Nl,sinal)

%%% Filtro de erro de predi��o (FEP)
%%% Estrutura: Linear
%%% Crit�rio: minimiza��o de L1/L2
%%% Busca: Particle Swarm
%%% entrada=vetor de entrada do FEP (deve ser um vetor linha)
%%% N=n�mero de coeficiente do filtro
%%% Np=n�mero de part�culas por vizinhan�a
%%% Nl=n�mero de vizinhan�as
%%% sinal=1 se o sinal a ser recuperado for super gaussiano, e -1 se o sinal for sub-gaussiano 

%%% entrada=vetor de entrada ( 1 coluna x Ns linhas) ou matriz ( Ne colunas
%%% x Ns linhas );

%%% Ne - corresponde ao n�mero de vetores que se deseja desconvoluir
%%% Ns - n�mero de amostras dos vetores que se deseja desconvoluir

%%% Autor: Kenji Nose Filho
%%% �ltima modifica��o: 23/06/2012

%% INICIALIZA��O

v_max=1; % velocidade m�xima
h=randn(Nl*Np,N); % iniciliza��o das part�culas
% h=h./(max(abs(h),[],2)*ones(N,1)'); % normaliza��o das part�culas pelo
% valor m�ximo
h=h ./ (sqrt(sum(h.^2,2)) * ones(1,N)); % normaliza��o das part�culas pela pot�ncia
v=zeros(Nl*Np,N); % inicializa��o das velocidades
melhor_custo_p=ones(Nl*Np,1)*1e10; % inicializa��o de pbest, melhor custo das particula
melhor_custo_l=ones(Nl,1)*1e10; % inicializa��o de lbest, melhor custo da vizinhan�a que � igual ao melhor custo global se Nl=1;
melhor_custo_g=1e10; % inicializa��o de gbest, melhor custo global

%%

for ii=1:N_iter
    f(ii)=1;
    ii/N_iter
    
    ll=1;
    for jj=1:Nl
        for kk=1:Np
            %% C�LCULO DA FUN��O CUSTO
            
            residuo_temp1=filter(h(ll,:),1,entrada); % c�lculo do res�duo de cada particula
            
            residuo_temp2=sum(abs(residuo_temp1));
            residuo_temp3=sqrt(sum(residuo_temp1.^2));
            
            residuo_temp4=residuo_temp2 ./ residuo_temp3;
                               
            custo(ll)=sinal * sum(residuo_temp4); % c�lculo do custo de cada particula
            custo_temp(kk)=custo(ll);
            h_temp(kk,:)=h(ll,:);
            
            %% ATUALIZA��O DE pbest
            
            if custo(ll)<=melhor_custo_p(ll)
                melhor_custo_p(ll)=custo(ll);
                melhor_h_p(ll,:)=h(ll,:);
            end
            
            ll=ll+1;
        end
        
        %% C�LCULO DO MELHOR CUSTO DA VIZINHAN�A jj NO INSTANTE ii
        indice_l=find(custo_temp==min(custo_temp),1);
        custo_l=custo_temp(indice_l);
        h_l=h_temp(indice_l,:);
        %% ATUALIZA��O DE lbest
        
        if custo_l<=melhor_custo_l(jj)
            melhor_custo_l(jj)=custo_l;
            melhor_h_l(jj,:)=h_l;
        end
        
    end
    
    %% C�LCULO DO MELHOR CUSTO GLOBAL NO INSTANTE ii
    indice_g=find(custo==min(custo),1);
    custo_g=custo(indice_g);
    h_g=h(indice_g,:);
    
    %% ATUALIZA��O DE gbest
    
    if custo_g<=melhor_custo_g
        melhor_custo_g=custo_g;
        melhor_h_g=h_g;
    end
    %%
    
    %% ATUALIZA��O DE v
    
    melhor_h_l_temp=[];
    for jj=1:Nl
        melhor_h_l_temp=[melhor_h_l_temp; ones(Np,1)*melhor_h_l(jj,:)];
    end
       
    atualiza_p=melhor_h_p-h;
    atualiza_l=melhor_h_l_temp-h;
    
%     v=v+2*rand(Np*Nl,N).*atualiza_p+2*rand(Np*Nl,N).*atualiza_l;
    v=v+2*atualiza_p+2*atualiza_l;
    
    v=v ./ (sqrt(sum(v.^2,2)) * ones(1,N));
    
%     for jj=1:Np
%         for kk=1:N
%             if abs(v(jj,kk))>v_max;
%                 v(jj,kk)=v_max*sign(v(jj,kk));
%             end
%         end
%     end
    
    erro_p(ii,:)=melhor_custo_p;
    erro_l(ii,:)=melhor_custo_l;
    erro_g(ii,:)=melhor_custo_g;
   
    if ii>1
        if sum(abs(erro_l(ii,:)-erro_l(ii-1,:)))<1e-2
            h=randn(Nl*Np,N);
            v=randn(Nl*Np,N);
            f(ii)=2;
        end         
    end
    
    h=h+v;
%     h=h./(max(abs(h),[],2)*ones(N,1)');
    h=h ./ (sqrt(sum(h.^2,2)) * ones(1,N)); % normaliza��o das part�culas pela pot�ncia
    
end

figure
subplot(131)
plot(erro_p)
subplot(132)
plot(erro_l)
subplot(133)
plot(erro_g)

% figure
% plot(f)

% melhor_h_g=melhor_h_g./melhor_h_g(1);
melhor_h_g = melhor_h_g ./ sqrt(sum(melhor_h_g.^2));
ii=find(melhor_h_g==max(abs(melhor_h_g)));
    
e=filter(melhor_h_g,1,entrada);

end



