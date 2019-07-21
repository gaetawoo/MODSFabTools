function map = zSynthesize(zCoef,maxDegree,rMin,rMax,maskVal)
%
% map = zSynthesize(zCoef,maxDegree,rMin,rMax,maskVal)
%
% Synthesizes 2-D map from Zernike polynomials.
%
% INPUT
% zCoef     Zernike coefficients (standard form, column vector)
% maxDegree of Zernike polynomials
% rMin      mask out area inside this radius
% rMax      mask out area outside this radius; normalizing radius
% maskVal   value to display for region outside aperture

% Set up arrays for map.
%
rows = 200;
cols = 200;
points = rows*cols;
dX = 2*rMax/cols;
dY = 2*rMax/rows;
Xvec = -rMax+dX/2:dX:rMax-dX/2; 
Yvec = rMax-dY/2:-dY:-rMax+dY/2; 
[X,Y] = meshgrid(Xvec,Yvec);
R = sqrt(X.*X+Y.*Y);
onMirror = R<=rMax & R>=rMin;
%
% Evaluate polynomials over map.
%
polys = size(zCoef);
if polys ~= (maxDegree+1)*(maxDegree+2)/2
    error('Wrong number of coefficients.')
end
rho = reshape(R/rMax,points,1);
theta = reshape(atan2(Y,X),points,1);
[zMatrix,nVec,elVec] = zEval(rho,theta,maxDegree);
map = reshape(zMatrix*zCoef,rows,cols);
map = onMirror .* map;
%
% Assign values to masked areas.
%
for i=1:rows
    for j=1:cols
        if ~onMirror(i,j)
            map(i,j) = maskVal;
        end
    end
end
