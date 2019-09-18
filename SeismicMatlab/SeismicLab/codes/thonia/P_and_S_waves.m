function P_and_S_waves

clear all
close all

f = 2;
T = 1/f
omega = 2*pi*f
vel = 3000

%xg =[-1.5,1.5]
%zg =[0,0]

% Set up array of depths
x1 = [0:150:5000]
x2 = x1
nx = length(x1)
k = 2*pi*f/vel

%Start time step loop
n=40;
M=moviein(n);
for it=1:n
  t=T*(it)/n;
  for ix=1:nx
     A1(ix) = 0.2*cos(k*x1(ix)-omega*t);
     A2(ix) = 0.0;
     x2(ix) = x1(ix)+100*cos(k*x1(ix)-omega*t);
  end
  
% Plot S-wave 
  plot(x1*0.001,A1+1.0,'o');
  
  hold on
  plot(x1*0.001,A1+1.3,'o');
  plot(x1*0.001,A1+1.15,'o');
  plot(x1*0.001,A1+0.85,'o');
  plot(x1*0.001,A1+0.7,'o');
  
  plot(x1(20)*0.001,A1(20)+1.30,'ro')
  plot(x1(20)*0.001,A1(20)+1.15,'ro')
  plot(x1(20)*0.001,A1(20)+1,'ro')
  plot(x1(20)*0.001,A1(20)+0.85,'ro')
  plot(x1(20)*0.001,A1(20)+0.70,'ro')
  
  plot(x1(1)*0.001,A1(1)+1.30,'ro')
  plot(x1(1)*0.001,A1(1)+1.15,'ro')
  plot(x1(1)*0.001,A1(1)+1,'ro')
  plot(x1(1)*0.001,A1(1)+0.85,'ro')
  plot(x1(1)*0.001,A1(1)+0.70,'ro')
  
% Plot P-wave
  plot(x2*0.001,A2-1.3,'o');
  plot(x2*0.001,A2-1.15,'o');
  plot(x2*0.001,A2-1.0,'o');
  plot(x2*0.001,A2-0.85,'o');
  plot(x2*0.001,A2-0.7,'o');
  
  plot(x2(20)*0.001,A2(20)-1.30,'ro');
  plot(x2(20)*0.001,A2(20)-1.15,'ro');
  plot(x2(20)*0.001,A2(20)-1,'ro');
  plot(x2(20)*0.001,A2(20)-0.85,'ro');
  plot(x2(20)*0.001,A2(20)-0.70,'ro');
  
  plot(x2(1)*0.001,A2(1)-1.30,'ro');
  plot(x2(1)*0.001,A2(1)-1.15,'ro');
  plot(x2(1)*0.001,A2(1)-1,'ro');
  plot(x2(1)*0.001,A2(1)-0.85,'ro');
  plot(x2(1)*0.001,A2(1)-0.7,'ro');
  
  hold off
  axis([-0.5,5,-2.2,2.2])
  title('P-waves and S-waves')
  ylabel ('A')
  xlabel ('(km)')
  
  M(:,it)=getframe;
end
  
movie(M,10)
  
  
print -dpsc  P-wave-and-S-wave.ps