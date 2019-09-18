function waves_attenuation


% Plots two waves propagating with absorbtion defined
% as 0.3 per wavelength

clear all
close all

f1 = 1
f2 = 4
omega1 = 2*pi*f1
omega2 = 2*pi*f2
vel = 3000

alpha = 0.3
lambda1 = vel/f1
lambda2 = vel/f2

%xg =[-1.5,1.5]
%zg =[0,0]

% Set up array of depths
ix=0
for xdum =0:50:5000;
   ix=ix+1;
   x(ix) = xdum;
end
nx = length(x)

k1 = 2*pi*f1/vel
k2 = 2*pi*f2/vel

%Start time step loop
n=50;
M=moviein(n);
for it=1:n
  t=(it)/n;
  for ix=1:nx
     A1(ix) =exp(-alpha*x(ix)/lambda1)*0.5*cos(k1*x(ix)-omega1*t);
     A2(ix) =exp(-alpha*x(ix)/lambda2)*0.5*cos(k2*x(ix)-omega2*t);
  end
  plot(x*0.001,A1+1,'b');
  axis([-0.5,5,-2.2,2.2])
  hold on
  plot(x(1)*0.001,A1(1)+1,'*')
  plot(x*0.001,A2-1,'b');
  plot(x(1)*0.001,A2(1)-1,'*')
  hold off
  title('f=1 and 4 Hz :    Absorbtion coefficient = 0.3 ')
  ylabel ('A')
  xlabel ('(km)')
  
  M(:,it)=getframe;
end
  
  movie(M,100)
  
  
print -dpsc  waves_1_and_4_hz_attenuation.ps