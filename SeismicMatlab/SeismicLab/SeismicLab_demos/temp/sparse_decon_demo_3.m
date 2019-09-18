%SPARSE_DECON_DEMO: l1 sparse-spike deconvolution
%
%
%            -----------------------------------
%              A demo for Sparse Deconvolution 
%
%                  M.D.Sacchi, SeismicLab
%            -----------------------------------
%
% The l1 regularization is used to estimate a sparse reflectivity 
% sequence. The problem involves minimizing a cost function
% J = ||Wr-d||^2 + mu * l1(r) where the first time is the misfit,
% convolution of wavelet W with reflectivity r should reproduce
% the data d, the second term is the regularization term used
% to impose sparsity on the solution.
%
% This demo uses the following SeismicLab functions:
%
%               readsegy.m
%               sparse_decon.m     
%               wigb.m
%

clear;
clear;
clc
close all;

setpath_seismiclab

% Read wavelet

% [w,H] = readsegy('wavelet_for_small_stack.su');
[w,H] = readsegy('min_phase_wavelet.su');
w=w';
w1=w;
kk=1;
for ii=length(w1):-1:1
    w2(kk)=w1(ii);
    kk=kk+1;
end
w=conv(w1,w2);    

NUM=w;
DEN=1;
SYS=tf(NUM,DEN,1,'variable','z^-1');
P=pole(SYS); Z=zero(SYS);
[H,T]=impulse(SYS);

figure,
subplot(121), zplane(Z,P)
subplot(122), plot(T,H)

% Read traces

% load data
% r=data;

r=full(sprandn(400,1,0.5))*ones(1,10);
r(401:500,:)=0;

for ii=1:size(r,2)
    d(:,ii)=filter(w,1,r(:,ii));
end

d=awgn(d,6);

N=25;

h0=zeros(1,N);
h0(ceil(N/2))=1;

[e1,h1]=fep_wiener(d,N,0,0);
[e2,h2]=fep_wiggins(d,h0,300,0.6,1,0);

[best1] = bestwiener(h1,N);
[best2] = bestwiener(h2,N);

figure, plot(best1./sqrt(sum(best1.^2))), hold on
plot(best2./sqrt(sum(best2.^2)),'r'), hold off

% Estimate the reflectivity using l1 (sparse) regularization

%  for ii=1:size(e2,1)
%      for jj=1:size(e2,2)
%          if abs(e2(ii,jj))>4e-6
%              e2(ii,jj)=sign(e2(ii,jj))*4e-6;
%          end
%      end
%  end

lag=44;
e2_temp=zeros(size(e2));
e2_temp(1:end-lag,:)=e2(lag+1:end,:);
 
figure,
subplot(141); wigb(d);
title('Traço')

subplot(142); wigb(e1);
title('Wiener')

subplot(143); wigb(e2_temp);
title('Wiggins')

subplot(144); wigb(r);
title('Refletividade')

max_iter = 100;
mu = 100;

[r0,dp0] = sparse_decon(d,w',mu,max_iter);
[r1,dp1] = sparse_decon(d,best1',mu,max_iter);
[r2,dp2] = sparse_decon(d,best2',mu,max_iter);

r0(end-20:end,:)=0;
r1(end-20:end,:)=0;
r2(end-20:end,:)=0;

figure,
subplot(141); wigb(r);
title('Refletividade')

subplot(142); wigb(r1);
title('Wiener')

subplot(143); wigb(r2);
title('Wiggins')

subplot(144); wigb(r0);
title('Refletividade Estimada')
