% generate mesh
% 1 left - 2 down - 3 right - 4 up
dx = (xMax - xMin) / ncellx;
dr = (rMaxtotal - rMin) / ncellrtotal;
[x,r] = meshgrid(linspace(xMin,xMax,ncellx+1),linspace(rMin,rMaxtotal,ncellrtotal+1));

for i = 1:ncellrtotal
    for j = 1:ncellx
        volumetotal(i,j) =  dr * dx; % 2-dim cylindrical volume
        faceAreatotal(i,j,1) =  dr;
        faceAreatotal(i,j,3) =  dr;
        faceAreatotal(i,j,4) =  dx;
        faceAreatotal(i,j,2) = dx;

    end
end


% normal vecter, colume vector
normal = [[-1;0],[0;-1],[1;0],[0;1]];

% neighbour index offset and distance
ioffset = [0,-1,0,1];
joffset = [-1,0,1,0];
d = [dx,dr,dx,dr];




disp('meshtotal has been generated*********')