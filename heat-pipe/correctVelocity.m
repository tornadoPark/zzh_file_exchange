for i=1:ncellr
    for j=1:ncellx
        u(i,j,:)= u(i,j,:)-VbyA(i,j)*(gradP(i,j,:)-gradPOld(i,j,:));
    end
end