dt = 2e-3;
fc = 40;
tw = -0.02:dt:0.03;
h  = (1-2*pi^2*fc^2*tw.^2).*exp(-pi^2*fc^2*tw.^2);
hh = hilbert(h);
th = -50*pi/180;
%entender essa filtragem
h  = cos(th)*real(hh)+sin(th)*imag(hh);
h  = h';