function [idx3] = select_traces(X,p_max)

[Ns,Nt] = size(X);
for ii=1:Nt
    X(:,ii)       = X(:,ii) - mean(X(:,ii));
end
norm_X  = sqrt(sum(X.*X));
Index_matrix = zeros(Nt);

for trace = 1:Nt
    Y       = X;
    norm_Y  = norm_X;
    x       = Y(:,trace)'*Y/norm_Y(trace); %Produto interno de um traço normalizado com todos os outros
    angulos = x./norm_Y;                   %Cos dos ângulos entre um traço e os outros (coef. correlação)
    [min_ang,ind] = min(abs(angulos));
    index_count   = 1;
    Index_matrix(trace,index_count) = ind;
    Y             = Y - Y(:,trace)*(x./norm_Y(trace)); %Ortonoliza em relação ao traço
    Y(:,trace)    = NaN;                               %Traço atual não conta mais
    while (min_ang < p_max)
        norm_Y   = sqrt(sum(Y.*Y));
        x        = Y(:,ind)'*Y/norm_Y(ind);        %Produto interno do traço eliminado com todos os outros
        angulos  = x./norm_Y;                      %Ângulos entre traço eliminado e os outros (coef. correlação)
        Y        = Y - Y(:,ind)*(x./norm_Y(ind));  %Ortonoliza em relação ao traço eliminado
        Y(:,ind) = NaN;
        [min_ang,ind] = min(abs(angulos));         %Escolhe o traço mais ortogonal
        index_count   = index_count + 1;
        Index_matrix(trace,index_count) = ind;
    end
end

ll = 1;
for ii = 1:Nt
    for jj = 1:Nt
        if Index_matrix(ii,jj) ~= 0
            idx0(ll,1) = max(ii,Index_matrix(ii,jj));
            idx0(ll,2) = -min(ii,Index_matrix(ii,jj));
            ll = ll+1;
        end
    end
end

idx0      = unique(idx0,'rows');
idx1      = idx0;
[~,i_idx] = sort(idx0(:,1));
idx1(:,1) = idx0(i_idx,1);
idx1(:,2) = idx0(i_idx,2);

idx2      = idx0;
[~,i_idx] = sort(idx1(:,2),'descend');
idx2(:,1) = idx1(i_idx,1);
idx2(:,2) = idx1(i_idx,2);

ll = 1;
for kk = 1:size(idx2,1)
    idx3(ll,:)   = [kk  idx2(kk,1)  idx2(kk,1) abs(idx2(kk,2)) abs(idx2(kk,2))];
    idx3(ll+1,:) = [kk  idx2(kk,2)  idx2(kk,2) abs(idx2(kk,1)) abs(idx2(kk,1))];
    idx3(ll+2,:) = [kk  idx2(kk,1)  idx2(kk,2) abs(idx2(kk,2)) abs(idx2(kk,1))];
    idx3(ll+3,:) = [kk  idx2(kk,2)  idx2(kk,1) abs(idx2(kk,1)) abs(idx2(kk,2))];
    ll           = ll+4;
end



