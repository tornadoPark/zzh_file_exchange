         %% T-equation
        AT = sparse(zeros(ncellx*ncellrtotal,ncellx*ncellrtotal));
        bT = zeros(ncellx*ncellrtotal,1);

        

        for i = 1:ncellrtotal
            for j = 1:ncellx
                imat = (i-1)*ncellx + j;
                for k = 1:4
                    jmat = (i+ ioffset(k) - 1)*ncellx + j+joffset(k);
                    if bTypetotal(i,j,k) == -1
                        uf(1) = 0.5* (utotal(i,j,1) + utotal(i+ioffset(k),j+joffset(k),1));
                        uf(2) = 0.5* (utotal(i,j,2) + utotal(i+ioffset(k),j+joffset(k),2));
                        F(i,j,k) = uf * normal(:,k) *faceAreatotal(i,j,k);

                        AT(imat,imat) = AT(imat,imat) + F(i,j,k)/2 +abs(F(i,j,k)/2) + alpha_f(i,j,k)*faceAreatotal(i,j,k)/d(k);
                        AT(imat,jmat) = F(i,j,k)/2 - abs(F(i,j,k)/2) - alpha_f(i,j,k)*faceAreatotal(i,j,k)/d(k);
                    elseif bTypetotal(i,j,k) == 0 && (k == 1 || k == 3)



                    elseif bTypetotal(i,j,k) == 0 && k == 2 && j<= ncellx*0.3
                        bT(imat)  = bT(imat) + alpha_f(i,j,k)*faceAreatotal(i,j,k)*norm(TGradBot);
                        
                        
                    elseif bTypetotal(i,j,k) == 0 && k == 2 &&  ncellx*0.3<j<=ncellx*0.5
                       
                        
                    elseif bTypetotal(i,j,k) == 0 && k == 2 && j> ncellx*0.5
                            AT(imat,imat) = AT(imat,imat) + alpha_f(i,j,k)*faceAreatotal(i,j,k)/(d(k) + Rhocp_sh*alpha_f(i,j,k)/h);
                            bT(imat) = bT(imat) + T_inf*alpha_f(i,j,k)*faceAreatotal(i,j,k)/(d(k) + Rhocp_sh*alpha_f(i,j,k)/h);   
                        
                    elseif bTypetotal(i,j,k) == 41 || bTypetotal(i,j,k) == 42 || bTypetotal(i,j,k) == 43
                       sigma_eff=0.3;
                        F(i,j,:) = sigma_eff * uLid *normal(:,k)*faceAreatotal(i,j,k);
                        AT(imat,imat) = AT(imat,imat) + F(i,j,k)/2 +abs(F(i,j,k)/2) + alpha_f(i,j,k)*faceAreatotal(i,j,k)/d(k);
                        bT(imat) = bT(imat) - (F(i,j,k)/2 - abs(F(i,j,k)/2) - alpha_f(i,j,k)*faceAreatotal(i,j,k)/d(k)) * TTop;%
                        
                    else
                        error('Error occurred.\n Undefined boundary!');
                    end
                end
                AT(imat,imat) = AT(imat,imat) + volumetotal(i,j) / dt;
                bT(imat) = bT(imat) + volumetotal(i,j)*TLast(i,j) /dt;

                % relax
                AT(imat,imat) = AT(imat,imat) / relaxT;
                bT(imat) = bT(imat) + (1-relaxT)*AT(imat,imat)*T(i,j);
            end
%               T(ncellrtotal-ncellr+1,:)=T(ncellrtotal-ncellr,:);
        end
        [TSol,TFlag] = bicgstab(AT,bT);
        T = (reshape(TSol,ncellx,ncellrtotal))';
      