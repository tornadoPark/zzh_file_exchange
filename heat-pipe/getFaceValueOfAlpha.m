function alpha_f = getFaceValueOfAlpha(T,alpha,ioffset,joffset,bTypetotal,ncellr,ncellrtotal,ncellx,epsilon)
% using harmonic interpolation
[ncellrtotal,ncellx] = size(alpha);
alpha_f = zeros(ncellrtotal,ncellx,4);

for i = 1:ncellrtotal
    for j = 1:ncellx
        for k = 1:4
            if bTypetotal(i,j,k) == -1
                alpha_f(i,j,k) = 2*alpha(i,j)*alpha(i+ioffset(k),j+joffset(k))...
                    /(alpha(i,j)+alpha(i+ioffset(k),j+joffset(k)));
            elseif bTypetotal(i,j,k) ~= -1
                alpha_f(i,j,k) = alpha(i,j);
            end
        end
    end
end

% for i = 1:ncellrtotal
%     for j = 1:ncellx
%         for k = 1:4
%             if bTypetotal(i,j,k) == -1
%                 alpha_f(i,j,k) = (alpha(i,j)+alpha(i+ioffset(k),j+joffset(k)))/2;
% 
%             elseif bTypetotal(i,j,k) ~= -1
%                 alpha_f(i,j,k) = alpha(i,j);
%             end
%         end
%     end
% end
% 
% T=T(ncellrtotal-ncellr+1:ncellrtotal,:);
% lambda_l = (124.67-(0.11381*T)+(5.5226*10e-5*T.^2)-(1.1842*10e-8*T.^3)); %l为液相
% lambda_s = ones(ncellr,ncellx)*16.3;  %s为wick区域固体
% lambda_sh= ones(ncellrtotal-ncellr,ncellx)*16.3;  %sh为壳体
% 
% cp_v = 200;
% cp_l = 1500+(3.432e-4.*T.^2)-(0.557.*T);
% cp_s = 1260;
% 
% rho_l = (219+(275.32*(1-T/2503.7))+(511.58*(1-T/2503.7).^0.5));
% rho_s = ones(ncellr,ncellx)*8000;
% rho_sh = 8000;
% 
% rhocp_l = rho_l.*cp_l ;
% rhocp_s = rho_s.*cp_s;
% rhocp_sh = zeros(ncellrtotal-ncellr,ncellx)+rho_sh *1260;
% 
% lambda_eff=(lambda_l.*((lambda_l+lambda_s)-(1-epsilon)*(lambda_l-lambda_s)))./((lambda_l+lambda_s)+(1-epsilon).*(lambda_l-lambda_s));
% rhocp_eff= (epsilon*rhocp_l+(1-epsilon)*rhocp_s);
% 
% lambda_ave = lambda_eff(1,:).*lambda_sh(10,:)./(lambda_eff(1,:) + lambda_sh(10,:));
% 
% for j = 1:ncellx
%     alpha_f(ncellrtotal - ncellr,j,4) =  lambda_ave(j)/rhocp_eff(1,j);
%     alpha_f(ncellrtotal - ncellr+1,j,2) =  lambda_ave(j)./rhocp_eff(1,j);
% end
end