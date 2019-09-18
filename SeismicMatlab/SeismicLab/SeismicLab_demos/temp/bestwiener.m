% Best Wiener

function [best, all, residue, res_wien, position,residual, H, R, corrs] = bestwiener(canal,ordem)

N = length(canal);
M = ordem;

NN = M + N - 1;


H = zeros(ordem, NN);

for pp = 1:ordem
   
   H(pp,pp:pp+N-1) = canal;
   
end

R = H*H';

for kk = 1:NN
   
   ZF = zeros(NN,1);
   ZF(kk) = 1;
     
   wiener(kk,:) = (inv(H*H')*H*ZF)';
   
   conjunta = H'*(wiener(kk,:))';
   
   residual (kk) = norm(conjunta - ZF);
   
   res_wien(kk) = 1 - wiener(kk,:)*H*ZF;
   
end

[residue,position] = min(residual);

best = wiener(position,:);
all = wiener;

corrs = H'*(inv(R))^2*H;

return

   

