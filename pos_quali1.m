clear all;
close all;
clc;
%criação do dado
dt = 4./1000;
%tmax = 4;
tmax = 3.996;
amp = [1, 0.8, 0.64, 0.512, 0.41]; %amplitudes dadas pelo professor Renato
f0 = 40;
L=1;
h = [0:50:450]; %10 traços, de 50 em 50
%profundidade dos eventos, em amostras: 250, 300, 450, 600, 800
%em milisegundos:
tau=[996,1196,1796,2396,3196];
%para segundos...
tau=tau.*(10^(-3));
%isto eu vou variando depois
%snr2=[128,64,32,24,18,15,12,9];
snr = 128 ;
v=[1500, 2000, 2250, 2500, 2700];
%%
%  OUT  d:         Data that consist of a superposition of reflections
%                  with hyerbolic  moveout (no avo)
%       t,h:       time and offset axes 
Nw             = length(ricker(f0,dt));
[d,h,t]=hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
%imagesc(h,t,d);
%wigb(d);
[Ns, Nx]= size(d);
L              = Ns-Nw+1;
A  = criaA(d,L,Nx); %passou no teste de dimensionalidade!
M=A'*A;
%DECOMPOSIÇÃO EM AUTOVETORES E AUTOVALORES DE M
[v_M,d_M]=eig(M);
%%%%%
%%
%teste para achar a função de refletividade a partir do hyperbolic events
[d2,h2,t2]=hyperbolic_events(dt,10^10,tmax,h,tau,v,amp,snr,L);
%wigb(d2)






