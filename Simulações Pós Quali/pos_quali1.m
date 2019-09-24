clear all;
close all;
clc;
setpath_seismiclab_win
%criação do dado
dt = 4./1000;
tmax = 3.996;
amp = [1, 0.8, 0.64, 0.512, 0.41]; %amplitudes dadas pelo professor Renato
f0 = 40;
L=1;
h = [0:200:1800]; %10 traços, de 100 em 100
%profundidade dos eventos, em amostras: 250, 300, 450, 600, 800
%em milisegundos:
tau=[996,1196,1796,2396,3196];
%para segundos...
tau=tau.*(10^(-3));
v=[1500, 1800, 2050, 2300, 2500];
[dr,hr,tr]=reflectivity_events(dt,tmax,h,tau,v,amp,L);%r pq é o real
%snr=[128,64,32,24,18,15,12,9]; linear, né?
%SNR eu vou variando depois
d3=conv2(dr,ricker(f0,dt));
d3= d3./std(d3(:)); %ESTÁ CERTO?
%RUÍDO
n              = randn(size(d3));
n              = n./std(n(:));
    
    SNR=128;
    %SNR_dB         = 10*log10(SNR);
    
    alpha          = 1./sqrt(SNR);
    n              = alpha*n;
d3_awgn         = d3+n;
Nw= length(ricker(f0,dt));
[Ns, Nx]= size(d3);
L              = Ns-Nw+1;
A  = criaA(d3_awgn,L,Nx); %passou no teste de dimensionalidade!
M=A'*A;
%DECOMPOSIÇÃO EM AUTOVETORES E AUTOVALORES DE M
[v_M,d_M]=eig(M);
 for ii =1:length(d_M)
        v_M(:,ii)   = v_M(:,ii)/norm(v_M(:,ii));
 end
%Cálculo de w
r         = dr(:)/norm(dr(:));
%%%%%%%%%%%%%%%%% DEVO NORMALIZAR AQUI MESMO?
w=v_M'*r;
%%
%%%%%%%%%%%%%%%%%%%%%% GRÁFICOS %%%%%%%%%%%%%%%%%%%%%%%%%%
% 1- com as amplitudes dos autovalores
% 2- com a norma 1 dos autovetores
% 3- com as amplitudes de w
% 4 - A refletividade original
% 5 - O traço sísmico com ruído
% 6 - O resultado obtido com o "melhor" autovetor
% 7 - O resultado obtido com o Wiener 
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% 1 amplitudes dos autovalores
x1=1:length(d_M);
figure(1), 
plot(x1,diag(d_M)), hold on
title('Magnitudes dos Autovalores')
xlim([x1(1) x1(end)]), grid on
xlabel('Autovalor')
ylabel('Magnitude')
% 2- com a norma 1 dos autovetores
%%
normas=zeros(1,length(v_M));
for vv=1:length(v_M);
    normas(1,vv)=norm(v_M(:,vv),1); 
end
figure(2), 
plot(x1,normas,'.'), hold on
title('Normas L1 dos Autovetores')
xlim([x1(1) x1(end)]), grid on
xlabel('Autovetor')
ylabel('Norma L1')
%%% 3- com as amplitudes de w
%% 
figure(3)
plot(x1,w)
hold on
title('Análise de Componentes do vetor w')
xlabel('Autovetor')
ylabel('Constante correspondente')
%% ATÉ AQUI
%%%%%%%%%%
%%%%%%%%%%
% 4 - A refletividade original e 5- O traço sísmico com ruído
%%
figure(4)
subplot(121), wigb(dr), title('Refletividades Originais')
subplot(122), wigb(d3), title('Traços com Ruído')
% 6 - O resultado obtido com o "melhor" autovetor E 
% 7 - O resultado obtido com o Wiener 
%%
%achar o melhor quando tiver o r 
for ii =1:length(d_M)
        %v_M(:,ii)   = v_M(:,ii)/norm(v_M(:,ii));
        rest      = v_M(:,ii);
        %r         = dr(:)/norm(dr(:)); %empilho o dr em vetor coluna e normalizo
        [rest,~] = delag1(r,rest,2*Nw);
        %rest      = rest(:);
        pcorr(ii) = rest'*r/(norm(rest)*norm(r));
    end
    
%     %flag
%     
% %     figure, plot(abs(pcorr))
    
     melhor= find(abs(pcorr)==max(abs(pcorr)));
     melhorp= max(abs(pcorr));
% 
figure(5)
plot(x1,pcorr,'.')
hold on
title('PCC de cada autovetor em relação à refletividade real')
xlabel('Autovetor')
ylabel('PCC')

%%
%VETOR ótimo:
w_otim=zeros(10000,1);
for hh=1:length(w)
    w_otim=w_otim+w(hh).*v_M(:,hh);
end
pcorr_w = w_otim'*r/(norm(w_otim)*norm(r));
figure(6)
y = [melhorp pcorr_w];
barh(y)
yticklabels({'Melhor Autovalor','W Calculado'})
xlabel('PCC')
ylabel('Metodo')







