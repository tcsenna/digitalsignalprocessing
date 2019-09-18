function [e_wiener,h]=fep_wiener(entrada_temp,N,D,E)

%%% Filtro de erro de predi��o (FEP)
%%% Estrutura: Linear
%%% Crit�rio: Wiener
%%% entrada=vetor de entrada do FEP (deve ser um vetor linha)
%%% N=n�mero de amostras do filtro de predi��o (coeficientes n�o nulos)
%%% D=n�mero de amostras do passo de predi��o (D-primeiros coeficintes nulos)
%%% E=% de branqueamento de R

entrada=entrada_temp;

h=zeros(1,N+D+1);

for jj=1:size(entrada,2)

    I=eye(1+N+D);

    kk=1;
    for ii=2+D:1+N+D
        U(:,kk)=filter(I(ii,:),1,entrada(:,jj));
        kk=kk+1;
    end

    R=U'*U;
    R=R+R(1,1)*E/100;

    h_temp=pinv(R)*U'*entrada(:,jj);

    h_wiener=[1 zeros(1,D) -h_temp'];
    
    h=h+h_wiener;
    
end

h=h./sqrt(sum(h.^2));

e_wiener=filter(h,1,entrada_temp);
end




