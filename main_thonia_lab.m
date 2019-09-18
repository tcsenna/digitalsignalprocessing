clc;
clear;
close all;
Dir1='C:\Users\Administrator\Desktop\codigos\SeismicMatlab\SeismicLab\codes\';
Dir2='C:\Users\Administrator\Desktop\codigos\SeismicMatlab\SeismicLab\';
Dir3='C:\Users\Administrator\Desktop\codigos\DeconvolutionToolbox\';
load dadocomp
%rcomp(end-30:end,:) = 0;
rtho=rcomp;
path(path, strcat(Dir1,'bp_filter'));
path(path, strcat(Dir1,'decon'));
path(path, strcat(Dir1,'fx'));
path(path, strcat(Dir1,'interpolation'));
path(path, strcat(Dir1,'kl_transform'));
path(path, strcat(Dir1,'radon_transforms'));
path(path, strcat(Dir1,'scaling_tapering'));
path(path, strcat(Dir1,'segy'));
path(path, strcat(Dir1,'seismic_plots'));
path(path, strcat(Dir1,'synthetics'));
path(path, strcat(Dir1,'velan_nmo'));
path(path, strcat(Dir1,'spectra'));

path(path, Dir2);
path(path, strcat(Dir2,'SeismicLab_demos'));
path(path, strcat(Dir2,'SeismicLab_data'));

path(path, Dir3);

%catches the default stream used by the random generator
% defaultStream = RandStream.getGlobalStream;
%
% %establishes the random seed
% stream_exp = RandStream('mt19937ar','seed',1); defaultStream.State = stream_exp.State;

trials = 1;
for kk = 1:trials;
    %r      = full(sprandn(100,3,0.2));
    %mais amostras
    %r      = full(sprandn(500,3,0.2));
    %Nx = 3;
    %r  = Data(1:250,1:Nx);
    
     SYS1=tf(rtho(:,1)',1);
     raiz1a=zero(SYS1);
     SYS2=tf(rtho(:,2)',1);
     raiz2a=zero(SYS2);
     SYS3=tf(rtho(:,3)',1);
     raiz3a=zero(SYS3);
%     %verificar sehá interseções entre os conjuntos de raízes
%     x1=intersect(raiz1a, raiz2a);
%     x2=intersect(raiz1a, raiz3a);
%     x3=intersect(raiz2a, raiz3a);
%     minLength = min([length(raiz1a), length(raiz2a), length(raiz3a)]);
%     % Removes any extra elements from the longer matrix
%     raiz1a = raiz1a(1:minLength);
%     raiz2a = raiz2a(1:minLength);
%     raiz3a = raiz3a(1:minLength);
    %problema de dimensionalidade que não sei como resolver
%     teste=poly(raiz1a);
%     r2=poly(raiz2a);
%     r3=poly(raiz3a);
%     r=[teste' r2' r3'];
    %r=[r1 r2 r3];
%     figure, imagesc(r);
    %number of samples in each r and number of reflectivity series used
    [Nsr,Nr] = size(rtho);
    %r  = [r; zeros(100,Nx)];
    dt = 2e-3;
    fc = 40;
    tw = -0.02:dt:0.03;
    h  = (1-2*pi^2*fc^2*tw.^2).*exp(-pi^2*fc^2*tw.^2);
    hh = hilbert(h);
    th = -50*pi/180;
    h  = cos(th)*real(hh)+sin(th)*imag(hh);
    h  = h';
    x  = [];
    
    
    for j=1:Nr
        x(:,j)=conv(h,rtho(:,j));
    end
    n              = randn(size(x));
    x              = x./std(x(:));
    n              = n./std(n(:));
    SNR            = 400000;
    SNR_dB         = 10*log10(SNR);
    alpha          = 1./sqrt(SNR);
    %x_awgn         = x+alpha*n;
    x_awgn         = x; % SNR = Infinito
    x_awgn         = x_awgn./std(x_awgn(:));
    Nw             = length(h);
    [Ns, Nx]       = size(x_awgn);
    L              = Ns-Nw+1;
    A              = criaA(x_awgn,L,Nx);
    M              = A'*A;
    [v,d]          = eig(M);
    %     figure(1)
    %     plot(diag(d))
    %rtho é o que eu quero achar, é meu alvo, a refletividade real.
    rtho         = rtho(:)/norm(rtho(:));
   
    for ii =1:length(d)
        v(:,ii)   = v(:,ii)/norm(v(:,ii));
        r_est     = v(:,ii);
        [r_est,~] = delag1(rtho,r_est,2*Nw);
        %r_est    = r_est(:);
        %pcorr(ii)= r_est(:)'*r(:)/(norm(r_est(:))*norm(r(:)));
        pcorr(ii) = r_est'*rtho/(norm(r_est)*norm(rtho));
    end   
       
%     figure, plot(abs(pcorr))
    
    melhor(kk)  = find(abs(pcorr)==max(abs(pcorr)));
    melhorp(kk) = max(abs(pcorr));
    
    %maior(kk)   = length(d);
    
    % v(:,melhor) = v(:,melhor)./norm(abs(v(:,melhor)),inf);
    % figure, plot(r(:,end)./norm(r(:,end),inf),'b','LineWidth',2), hold on
    % plot(sign(pcorr(melhor))*v(:,melhor),'r','LineWidth',2)
end
%figurasthonia
figure, hist(melhor)
xlabel('Autovalor')
ylabel('% de melhores soluções em relação a cada autovetor/autovalor')
figure, hist(melhorp)
xlabel('PCC')
ylabel('% de maior PCC por autovetor/autovalor')
figure, plot(melhor,melhorp,'o')
xlabel('Autovalor')
ylabel('PCC')
figure, 
subplot(131), plot(raiz1a,'o'), hold,  plot(raiz2a,'rx')
subplot(132), plot(raiz1a,'o'), hold,  plot(raiz3a,'g*')
subplot(133), plot(raiz2a,'rx'), hold,  plot(raiz3a,'g*')
%%começa o do kenji
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deconvolution with multichannel algorithm
%%
Nw     = length(h);
Niter  = 1000;     % 200 (total number of iterations of the algorithm)
Nline  = 5;        % 5 (number of iteration of the line search)
epi    = 0.0005;   % 0.0005

%% Deconvolution with multichannel algorithm (reduced-smbd)

pcorr2  = [0.001 0.01];    % 0.001
lambda = [2.0   2.0 ];    % 2.0;

for kk=1:length(pcorr2);

[r_kenji,~,custo_multicanal20,custo_multicanal21,custo_multicanal22,perc(1,kk)] = mult_decon_reduced(x_awgn,Niter,Nline,epi,lambda(1,kk),Nw,pcorr2(1,kk));
[r_kenji,r_multicanal2_lag] = delag1(rcomp,r_kenji,2*Nw);

figure(),
subplot(234), plot(custo_multicanal21), title('J1')
subplot(235), plot(custo_multicanal22), title('lambda*J2')
subplot(236), plot(custo_multicanal20), title('J1+lambda*J2')

figure,
subplot(131), wigb(x_awgn)
subplot(132), wigb(r_kenji)
subplot(133), wigb(rcomp)

figure,
plot(r_kenji(:,1),'b'), hold on
plot(rcomp(:,1),'r-x'), hold off

pcorr_rr2(1,kk) = corr(rcomp(:),r_kenji(:));
end

perc
pcorr_rr2