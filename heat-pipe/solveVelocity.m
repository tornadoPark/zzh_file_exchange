 %% u-equation

        A = sparse(zeros(ncellx*ncellr,ncellx*ncellr));
        b = zeros(ncellx*ncellr,2);

        for i = 1:ncellr
            for j = 1:ncellx
                imat = (i-1)*ncellx + j;
                for k = 1:4
                    jmat = (i+ioffset(k)-1)*ncellx + j+joffset(k); % determine which neighbour delt with in this cycle
                    if bType(i,j,k) == -1
                        uf(1) = 0.5* (u(i,j,1) + u(i+ioffset(k),j+joffset(k),1));
                        uf(2) = 0.5* (u(i,j,2) + u(i+ioffset(k),j+joffset(k),2));
                        F(i,j,k) = uf*normal(:,k)*faceArea(i,j,k);

                        A(imat,jmat) = F(i,j,k)/2 - abs(F(i,j,k)/2) -nu*faceArea(i,j,k)/d(k);
                        A(imat,imat) = A(imat,imat) + (F(i,j,k)/2 +abs(F(i,j,k)/2) + nu*faceArea(i,j,k)/d(k));

                    elseif bType(i,j,k) == 0 
                        F(i,j,k) = 0;
                        A(imat,imat) = A(imat,imat) + nu*faceArea(i,j,k)/(d(k)/2);

                    elseif bType(i,j,k) == 41 %   evap-outlet
                        uf(1) = u(i,j,1);
                        uf(2) = u(i,j,2);
                        F(i,j,k) = uf*normal(:,k)*faceArea(i,j,k);
                        A(imat,imat) = A(imat,imat) + (F(i,j,k)/2 +abs(F(i,j,k)/2));
                        
                    elseif bType(i,j,k) == 42
                         F(i,j,k) = 0;
                        A(imat,imat) = A(imat,imat) + nu*faceArea(i,j,k)/(d(k)/2);
                       
                        
                    elseif bType(i,j,k) == 43 %  cond-inlet
%                         uLid = getuLid(i,ncellr,dr,Q_evap,l_cond,rho_l,hfg_vapor);
                        F(i,j,k) = uLid*normal(:,k)*faceArea(i,j,k);
                        
                        A(imat,imat) = A(imat,imat) + (F(i,j,k)/2 +abs(F(i,j,k)/2) + nu*faceArea(i,j,k)/d(k));
                        b(imat,1) = b(imat,1) + F(i,j,k)*uLid(1) +nu*uLid(1)*faceArea(i,j,k)/(d(k)/2);
                        b(imat,2) = b(imat,2) + F(i,j,k)*uLid(2) +nu*uLid(2)*faceArea(i,j,k)/(d(k)/2);

                    else
                        error('Error occurred.\n Undefined boundary!');
                    end
                end

                
                b(imat,1) = b(imat,1) - gradP(i,j,1)*volume(i,j);
                b(imat,2) = b(imat,2) - gradP(i,j,2)*volume(i,j);

                %relax
                A(imat,imat) = (A(imat,imat) +(nu * epsilon / kporous)*volume(i,j)+volume(i,j)/dt)/relaxU;
                %A(imat,imat) = (A(imat,imat)+ volume(i,j)/dt)/relaxU; % till now, diagnol element has been changed(relaxed)!!!
                b(imat,1) =b(imat,1) +u(i,j,1)*A(imat,imat) * (1-relaxU) + volume(i,j)*uLast(i,j,1)/dt;
                b(imat,2) =b(imat,2) +u(i,j,2)*A(imat,imat) * (1-relaxU) + volume(i,j)*uLast(i,j,2)/dt;
            end
        end

        %    aP = reshape(diag(A),ncellx,ncellr)'; % diagnol element of coeff matrix A
        for i = 1:ncellr
            for j = 1:ncellx
                imat = (i-1)*ncellx + j;
                aP(i,j) = A(imat,imat);% aP is under-relax value, due to A(imat,imat) has been relaxed.
            end
        end

        VbyA = volume./(aP); % pesdo conductive coeff of equ of pressure = VP/aP = volvin = volume velocity inverse = V* rAU from Openfoam, has been relaxed due to aP
        %     ux = (reshape(A\b(:,1),ncellx,ncellr))';% cause the inner cycle index is j, so the reshape shuould be in order of ncellx(for j) and then ncellr.
        %     ur = (reshape(A\b(:,2),ncellx,ncellr))';
        [u1Sol, u1Flag] = bicgstab(A,b(:,1),1e-9,100);
        [u2Sol, u2Flag] = bicgstab(A,b(:,2),1e-9,100);
        ux = (reshape(u1Sol,ncellx,ncellr))';
        ur = (reshape(u2Sol,ncellx,ncellr))';

        u(:,:,1) = ux;
        u(:,:,2) = ur;


        %     disp('momentum equation is solved********')