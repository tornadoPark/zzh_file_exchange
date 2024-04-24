%% p-equation

Ap = sparse(zeros(ncellx*ncellr,ncellx*ncellr));
bp = zeros(ncellx*ncellr,1);
for i = 1:ncellr
    for j = 1:ncellx
        imat = (i-1)*ncellx+j;
        for k = 1:4
            jmat = (i-1+ioffset(k))*ncellx + j+joffset(k);
            if bType(i,j,k) == -1
                VbyA_f(i,j,k) = (VbyA(i,j) + VbyA(i+ioffset(k),j+joffset(k)))/2.0;% pesdo condc coeff of equ pressure interpolated to face
                
                uP = [u(i,j,1),u(i,j,2)]; % predict value
                uN = [u(i + ioffset(k) ,j + joffset(k),1),u(i + ioffset(k) ,j + joffset(k),2)];
                gradPP = [gradP(i,j,1),gradP(i,j,2)]; % pressure on cell centroid
                gradPN = [gradP(i + ioffset(k) ,j + joffset(k),1),gradP(i + ioffset(k) ,j + joffset(k),2)];
                VbyAP = VbyA(i,j);
                VbyAN = VbyA(i + ioffset(k) ,j + joffset(k));
                sn = normal(:,k);  % surface normal
                %
                HbyA_f(i,j,k) = ((uP + VbyAP*gradPP + uN + VbyAN*gradPN)/2.0)*sn;
                %                     HbyA_f(i,j,k) = (HbyA(i,j) + HbyA(i+ioffset(k),j+joffset(k))/2.0;
                
                Ap(imat,imat) = Ap(imat,imat) - VbyA_f(i,j,k)*faceArea(i,j,k)/d(k);
                Ap(imat,jmat) = VbyA_f(i,j,k)*faceArea(i,j,k)/d(k);
                
                bp(imat) = bp(imat) + HbyA_f(i,j,k)*faceArea(i,j,k);
                
            elseif bType(i,j,k) == 0 
                Ap(imat,imat) = Ap(imat,imat) - 0;
                bp(imat) = bp(imat)+0;
                % bp(imat) = bp(imat) + VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k);
                %  VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k)
                
            elseif bType(i,j,k) == 41 %   evap-outlet
                VbyA_f(i,j,k) = VbyA(i,j);
                uP = [u(i,j,1),u(i,j,2)];
                
                VbyAP = VbyA(i,j);
                gradPP = [gradP(i,j,1),gradP(i,j,2)]; % pressure on cell centroid
                
                sn = normal(:,k);
                HbyA_f(i,j,k) = (uP + VbyAP*gradPP)*sn;
                
                Ap(imat,imat) = Ap(imat,imat) - VbyA_f(i,j,k)*faceArea(i,j,k)*2/d(k);
                bp(imat) = bp(imat) + HbyA_f(i,j,k)*faceArea(i,j,k) - pTop*(VbyA_f(i,j,k)*faceArea(i,j,k)*2/d(k));
                
            elseif bType(i,j,k) == 42 %
                Ap(imat,imat) = Ap(imat,imat) - 0;
                bp(imat) = bp(imat)+0;
                
            elseif bType(i,j,k) == 43 %   cond-inlet
                VbyAP = VbyA(i,j);
                gradPP = [gradP(i,j,1),gradP(i,j,2)]; % pressure on cell centroid
                
                sn = normal(:,k);
%                 uLid = getuLid(i,ncellr,dr);
                HbyA_f(i,j,k) = (uLid + VbyAP*gradPP)*sn;

                Ap(imat,imat) = Ap(imat,imat) - 0;
                bp(imat) = bp(imat) + HbyA_f(i,j,k)*faceArea(i,j,k);
                % bp(imat) = bp(imat) + VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k);
                %  VbyA(i,j)*([gradP(i,j,1),gradP(i,j,2)]*normal(:,k))*faceArea(i,j,k)
                
                
            else
                error('Error occurred. \n Undefined boundary!');
            end
        end
    end
end

pOld = p;
%     p = (reshape(Ap\bp,ncellx,ncellr))';
[pSol, pFlag] = bicgstab(Ap,bp,1e-9,1000);
p = (reshape(pSol,ncellx,ncellr))';
pUnCor = p;
p = pOld + relaxP *(p-pOld); % Jasak thesis equ 3.145

%     disp('continuity equation is solved********')