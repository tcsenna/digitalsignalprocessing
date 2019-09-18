clear;
close all;

%catches the default stream used by the random generator
% defaultStream = RandStream.getGlobalStream;
%
% %establishes the random seed
% stream_exp = RandStream('mt19937ar','seed',1); defaultStream.State = stream_exp.State;
setpath_seismiclab_dspcom;
trials = 1;
for kk = 1:trials;
    r1      = full(sprandn(500,3,0.2));
    %     figure, imagesc(r);
    %number of samples in each r and number of reflectivity series used
    [Nsr,Nr] = size(r);
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
        x(:,j)=conv(h,r(:,j));
    end
    n              = randn(size(x));
    x              = x./std(x(:));
    n              = n./std(n(:));
    %%
    SNR            = 4000;
    SNR_dB         = 10*log10(SNR);
    %%
    alpha          = 1./sqrt(SNR);
    n              = alpha*n;
    x_awgn         = x;
    x_awgn         = x_awgn./std(x_awgn(:));
    Nw             = length(h);
    [Ns, Nx]       = size(x_awgn);
    L              = Ns-Nw+1;
    A  = criaA(x_awgn,L,Nx);
    M=A'*A+eps;
    [u,d,v,flag] = svds(M,length(M),'smallest','Tolerance',1e-6);
%     figure(1)
%     plot(diag(d))
    
    for ii =1:length(d)
        v(:,ii)   = v(:,ii)/norm(v(:,ii));
        rest      = v(:,ii);
        r         = r(:)/norm(r(:));
        rest      = rest(:);
        pcorr(ii) = rest(:)'*r(:)/(norm(rest(:))*norm(r(:)));
    end
    
    flag
    
%     figure, plot(abs(pcorr))
    
    melhor(kk)  = find(abs(pcorr)==max(abs(pcorr)))
    melhorp(kk) = max(abs(pcorr))
    maior(kk)   = length(d);
    
    % v(:,melhor) = v(:,melhor)./norm(abs(v(:,melhor)),inf);
    % figure, plot(r(:,end)./norm(r(:,end),inf),'b','LineWidth',2), hold on
    % plot(sign(pcorr(melhor))*v(:,melhor),'r','LineWidth',2)
end

figure, hist(melhor)
figure, hist(melhorp)

figure, plot(melhor,melhorp,'o')