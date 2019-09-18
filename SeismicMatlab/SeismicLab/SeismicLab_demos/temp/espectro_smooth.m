function [modulo_db]=espectro_smooth(data,fa)

%%% Funï¿½ï¿½o que calcula e plota o espectro de frequencia de um sinal data,
%%% com frequï¿½ncia de amostragem fa.

N=size(data,1);     % nï¿½mero de amostras
t=(0:N-1)/fa;       % vetor tempo
w=(0:N-1)*fa/(N-1); % vetor frequï¿½ncia

for ii=1:size(data,2)

    f(:,ii)=fft(data(:,ii));

%     modulo_db(:,ii)=20*log10(abs(f(:,ii))+eps);
    modulo_db(:,ii)=(abs(f(:,ii))+eps);
%     modulo_db(:,ii)=modulo_db(:,ii) / max(modulo_db(:,ii));
end
% % modulo_db=abs(f);
% temp=angle(f);
% fase=unwrap(angle(f));
% fase=fase*180/pi;
% for ii=1:length(fase)
%     if fase(ii)>180
%         temp1(ii)=fase(ii)/180;
%         temp2(ii)=fase(ii)-floor(temp1(ii))*180;
%     end
%     if fase(ii)<-180
%         temp1(ii)=fase(ii)/-180;
%         temp2(ii)=fase(ii)+floor(temp1(ii))*180;
%     end
% end
% fase=temp2;
% 
% N2=floor(N/2);
% 
% figure,
% subplot(211)
% plot(w(1:N2),modulo_db(1:N2))
% axis([w(1) w(N2) min(modulo_db(1:N2)) 1.1*max(modulo_db(1:N2))])
% ylabel('Mï¿½dulo [dB]')
% xlabel('Frequï¿½ncia [HZ]')
% subplot(212)
% plot(w(1:N2),fase(1:N2))
% axis([w(1) w(N2) min(fase(1:N2)) 1.1*max(fase(1:N2))])
% ylabel('Fase [rad]')
% xlabel('Frequï¿½ncia [HZ]')

media=20*log10(mean(modulo_db,2));
N2=floor(N/2);
figure,
plot(w(1:N2),(media(1:N2)))
axis([w(1) w(N2) min(media(1:N2)) 1.1*max(media(1:N2))])
ylabel('Módulo [dB]')
xlabel('Frequência [HZ]')

end