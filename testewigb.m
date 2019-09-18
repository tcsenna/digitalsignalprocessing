clc;
clear;
close all;
load dadozero_tho

setpath_seismiclab_dspcom;

%catches the default stream used by the random generator
defaultStream = RandStream.getGlobalStream;

%establishes the random seedb
stream_exp = RandStream('mt19937ar','seed',1); defaultStream.State = stream_exp.State;

Nx = 15;
r  = dadozero(1:250,1:Nx);
r  = [r; zeros(100,Nx)];
dt =2e-3;

fc = 40;
%eixo temporal da wavelet
tw = -0.02:dt:0.03;
%pura e simplesmente a wavelet
h  = (1-2*pi^2*fc^2*tw.^2).*exp(-pi^2*fc^2*tw.^2);
%pq pegou a transf hilbert da wavelet? nao entendi
hh = hilbert(h);
%tb n entendi isso de th
th = -50*pi/180;
h  = cos(th)*real(hh)+sin(th)*imag(hh);
h  = h';
%dei um help filter e não entendi o 1 aqui no meio.achei q era um vetor
%x é o dado r filtrado com o filtro h
x  = filter(h,1,r);

%% Aditive White Gaussian Noise
n              = randn(size(r)); %que é 350x15
%fez isso pq agora o desvio padrão é um
x              = x./std(x(:));
%fez o mesmo pra n
n              = n./std(n(:));
SNR            = 1000;
SNR_dB         = 10*log10(SNR);
alpha          = 1./sqrt(SNR); %deu 0.0316

n              = alpha*n; %pq ele fez isso?
x_awgn         = x + n;
[x_awgn,d_lag] = delag1(r,x_awgn,length(h));
x_awgn         = x_awgn./std(x_awgn(:));

%% Deconvolution with multichannel algorithm

Nw     = length(h);
Niter  = 1000;     % 200 (total number of iterations of the algorithm)
Nline  = 5;        % 5 (number of iteration of the line search)
epi    = 0.0005;   % 0.0005
lambda = 3.0;      % 3.0
tic
[r_multicanal1,custo_multicanal10,custo_multicanal11,custo_multicanal12] = mult_decon(x_awgn,Niter,Nline,epi,lambda,Nw);
time_total = toc
r_multicanal1                     = [zeros(25,Nx); r_multicanal1];
[r_multicanal1,r_multicanal1_lag] = delag1(r,r_multicanal1,2*Nw);

