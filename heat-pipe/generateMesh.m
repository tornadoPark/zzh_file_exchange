% generate mesh
% 1 left - 2 down - 3 right - 4 up
dx = (xMax - xMin) / ncellx;
dr = (rMax - rMin) / ncellr;
[x,r] = meshgrid(linspace(xMin,xMax,ncellx+1),linspace(rMin,rMax,ncellr+1));

for i = 1:ncellr
    for j = 1:ncellx
        volume(i,j) =  dr * dx; % 2-dim cylindrical volume
        faceArea(i,j,1) =  dr;
        faceArea(i,j,3) =  dr;
        faceArea(i,j,4) =  dx;
        faceArea(i,j,2) = dx;

    end
end


% normal vecter, colume vector
normal = [[-1;0],[0;-1],[1;0],[0;1]];

% neighbour index offset and distance
ioffset = [0,-1,0,1];
joffset = [-1,0,1,0];
d = [dx,dr,dx,dr];




disp('mesh has been generated*********')