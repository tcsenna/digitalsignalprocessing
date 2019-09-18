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

 clc;
 clear;
 close all;
 
 setpath_seismiclab
% Read traces

 [d,Hs]=readsegy('small_stack.su');

% Extract sampling interval and cdp number

  dtsec = Hs(1).dt/1000/1000;
  cdp = [Hs.cdp];

% Read wavelet

 [w,Hw] = readsegy('wavelet_for_small_stack.su');
 w=awgn(w,24);
 figure, plot(w)
% Estimate the reflectivity using l1 (sparse) regularization

 max_iter = 20;
 mu = 0.1;
 
 [best] = bestwiener(w,100);
 
 r = filter(best,1,d);
 
 figure, stem(conv(w,best))

%  [r,dp] = sparse_decon(d,w,mu,max_iter);
% 
% Display results

 figure,;

 [nt,ncdp] = size(d);
 t_axis =[0:1:nt-1]*dtsec;

 subplot(131); wigb(d,2,cdp,t_axis);
 xlabel('cdp');
 ylabel('Time (s)');

 subplot(132); wigb(r,2,cdp,t_axis);
 xlabel('cdp');

%  subplot(133); wigb(dp,2,cdp,t_axis);
%  xlabel('cdp');

aux_d = d(:);
aux_r = r(:);

figure, plot(aux_r,filter([0 1],1,aux_r),'x')

