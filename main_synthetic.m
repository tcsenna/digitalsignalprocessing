%%% Kenji Nose Filho - kenjinose@yahoo.com.br
%%% 05/04/2017

clc;
clear;
close all;

load synthetic2

setpath_seismiclab_dspcom;
%%
%catches the default stream used by the random generator
defaultStream = RandStream.getGlobalStream;

%establishes the random seedb
stream_exp = RandStream('mt19937ar','seed',1); defaultStream.State = stream_exp.State;

Nx = 15;
r  = Data(1:250,1:Nx);
r  = [r; zeros(100,Nx)];
dt = 2e-3;

fc = 40;
tw = -0.02:dt:0.03;
h  = (1-2*pi^2*fc^2*tw.^2).*exp(-pi^2*fc^2*tw.^2);
hh = hilbert(h);
%a parte real de hh é igual a h. o diferente eh o componente imaginario
th = -50*pi/180;
%entender essa filtragem
h  = cos(th)*real(hh)+sin(th)*imag(hh);
%h que tinha 1x26 passa a ter 26x1
h  = h';
x  = filter(h,1,r);
%agora x é uma parte do dado com zeros embaixo que tá filtrada
% Aditive White Gaussian Noise
n              = randn(size(r));
x              = x./std(x(:));
n              = n./std(n(:));
%assim o std de x e de n passam a ser 1
SNR            = 1000;
SNR_dB         = 10*log10(SNR);
alpha          = 1./sqrt(SNR);
n              = alpha*n;
%foi adicionado a x o ruido branco gaussiano
x_awgn         = x + n;
%ENTENDER ESSA PARTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[x_awgn,d_lag] = delag1(r,x_awgn,length(h));
x_awgn         = x_awgn./std(x_awgn(:));

%Deconvolution with multichannel algorithm
%O programa original do sacchi
Nw     = length(h);
Niter  = 1000;     % 200 (total number of iterations of the algorithm)
Nline  = 5;        % 5 (number of iteration of the line search)
epi    = 0.0005;   % 0.0005
lambda = 3.0;      % 3.0
%%
tic
[r_multicanal1,custo_multicanal10,custo_multicanal11,custo_multicanal12] = mult_decon(x_awgn,Niter,Nline,epi,lambda,Nw);
time_total = toc
r_multicanal1                     = [zeros(25,Nx); r_multicanal1];
[r_multicanal1,r_multicanal1_lag] = delag1(r,r_multicanal1,2*Nw);

%% Deconvolution with multichannel algorithm (reduced-smbd)

pcorr  = [0.001 0.005 0.02 0.05 0.1];        % [0.001 0.005 0.02 0.05 0.1];
lambda = [2.0   2.0   3.0  3.0  3.0];        % [2.0   2.0   3.0  3.0  3.0];

for kk=1:length(pcorr);

tic
[r_multicanal2,~,custo_multicanal20,custo_multicanal21,custo_multicanal22,perc(1,kk)] = mult_decon_reduced(x_awgn,Niter,Nline,epi,lambda(1,kk),Nw,pcorr(1,kk));
time_reduced(1,kk) = toc
r_multicanal2                     = [zeros(25,Nx); r_multicanal2];
[r_multicanal2,r_multicanal2_lag] = delag1(r,r_multicanal2,2*Nw);

figure(300+kk),
subplot(231), plot(custo_multicanal11), title('J1')
subplot(232), plot(custo_multicanal12), title('lambda*J2')
subplot(233), plot(custo_multicanal10), title('J1+lambda*J2')
subplot(234), plot(custo_multicanal21), title('J1')
subplot(235), plot(custo_multicanal22), title('lambda*J2')
subplot(236), plot(custo_multicanal20), title('J1+lambda*J2')

[P_r,f]        = smooth_spectrum(r,dt,0,'db');
P_dn           = smooth_spectrum(x_awgn,dt,0,'db');
P_multicanal1  = smooth_spectrum(r_multicanal1,dt,0,'db');
P_multicanal2  = smooth_spectrum(r_multicanal2,dt,0,'db');
perc(1,kk)
pcorr_rr1(1,kk)  = corr(r(:),r_multicanal1(:))
pcorr_rr2(1,kk)  = corr(r(:),r_multicanal2(:))
pcorr_r1r2(1,kk) = corr(r_multicanal1(:),r_multicanal2(:))

end

figure(1), 
semilogx(pcorr,pcorr_r1r2,'k','LineWidth',2), hold on
axis([pcorr(1) pcorr(end) 0.97 1]), grid on
xlabel('\rho_{max}')
ylabel('\rho_{smbd x reduced-smbd}')

figure(2),
semilogx(pcorr,pcorr_rr1,'k--','LineWidth',2), hold on
semilogx(pcorr,pcorr_rr2,'k','LineWidth',2), hold off
axis([pcorr(1) pcorr(end) 0.78 0.795]), grid on
xlabel('\rho_{max}')
legend('\rho_{real x smbd}','\rho_{real x reduced-smbd}')

figure(3),
semilogx(pcorr,perc,'k','LineWidth',2),
axis([pcorr(1) pcorr(end) 10 100]), grid on
xlabel('\rho_{max}')
ylabel('% of total combinations')

position1 = [10 51 300 300]
set(1,'Position',position1)
set(2,'Position',position1)
set(3,'Position',position1)

%%

% figure(300)
% subplot(231), plot(custo_multicanal11), title('J1')
% subplot(232), plot(custo_multicanal12), title('lambda*J2')
% subplot(233), plot(custo_multicanal10), title('J1+lambda*J2')
% 
% subplot(234), plot(custo_multicanal21), title('J1')
% subplot(235), plot(custo_multicanal22), title('lambda*J2')
% subplot(236), plot(custo_multicanal20), title('J1+lambda*J2')

figure(4), wigb(r,1,1:Nx,(0:349)*dt), ylabel('Time(s)'), xlabel('Trace number')
figure(5), wigb(x_awgn,1,1:Nx,(0:349)*dt), ylabel('Time(s)'), xlabel('Trace number')
figure(6), wigb(r_multicanal1,1,1:Nx,(0:349)*dt), ylabel('Time(s)'), xlabel('Trace number')
figure(7), wigb(r_multicanal2,1,1:Nx,(0:349)*dt), ylabel('Time(s)'), xlabel('Trace number')

position1 = [10 51 400 600];

set(4,'Position',position1)
set(5,'Position',position1)
set(6,'Position',position1)
set(7,'Position',position1)

% figure(7), plot(f,P_r,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% figure(8), plot(f,P_dn,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% figure(9), plot(f,P_huber,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% figure(10), plot(f,P_wiener,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% figure(11), plot(f,P_multicanal1,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% figure(12), plot(f,P_multicanal2,'k','Linewidth',2), axis([0 f(end) -20 0]); grid on, xlabel('Freq. (Hz)')
% 
% position2 = [10 51 400 150];
% 
% set(7,'Position',position2)
% set(8,'Position',position2)
% set(9,'Position',position2)
% set(10,'Position',position2)
% set(11,'Position',position2)
% set(12,'Position',position2)