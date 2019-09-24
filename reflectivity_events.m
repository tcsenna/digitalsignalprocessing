function [d,h,t] = reflectivity_events(dt,tmax,h,tau,v,amp,L);
%code to generate the original reflectivity
%[d2,h2,t2]=reflectivity_events(dt,tmax,h,tau,v,amp,L);
%  [d] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
%
%  IN   dt:        sampling interval in secs
%       tmax:      maximun time of the simulation in secs
%       h:         vector of offsets in meters
%       tau,v,amp: vectors of intercept, rms velocities
%                  and amplitude of each linear event
%                  (v is in m/s and tau in secs)
%       L:         The random noise is average over L samples
%                  to simulate band-pass noise (L=1 means no averaging)
%
%  OUT  d:         Data that consist of a superposition of reflections
%                  with hyerbolic  moveout (no avo)
%       t,h:       time and offset axes 
%
%  Example with default parameters:
%
%    [d,h,t] = hyperbolic_events; imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%

 nt = floor(tmax/dt)+1;
 nfft = 4*(2^nextpow2(nt));
 n_events = length(tau);
 nh = length(h);
 wavelet = 1;
 nw = length(wavelet);
 W = fft(wavelet,nfft);
 W = ones(nfft,1);
 D = zeros(nfft,nh);
 i = sqrt(-1);

% Important: the following lines is to have the maximum of the Ricker
% wavelet at the right intercept time

 delay = dt*(floor(nw/2)+1);

for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(  sqrt(tau(k)^2 + (h/v(k)).^2) - delay));
    D(ifreq,:) = D(ifreq,:) +amp(k)* W(ifreq)*Shift;
  end
end

% Apply w-domain symmetries

 for ifreq=2:nfft/2
  D(nfft+2-ifreq,:) = conj(D(ifreq,:));
 end 

 d = ifft(D,[],1);
 d = real(d(1:nt,:));
 if nargout>1;
  t = (0:1:nt-1)*dt;
 end;

 return;
end

