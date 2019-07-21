function zDisplay(zCoef,maxDegree,rMin,rMax,palette,scale,maskVal)
%
% zDisplay(zCoef,maxDegree,rMin,rMax,palette,scale,maskVal)
%
% Displays map made of Zernike polynomials.
%
% INPUT
% zCoef     Zernike coefficients (standard form, column vector)
% maxDegree of Zernike polynomials
% rMin      mask out area inside this radius
% rMax      mask out area outside this radius; normalizing radius
% palette   Matlab colormap (string)
% scale     [min,max] for color map 
% maskVal   value to display for region outside aperture

% Set up arrays for map.
%
rows = 200;
cols = 200;
points = rows*cols;
dX = 2*rMax/cols;
Xvec = -rMax+dX/2:dX:rMax-dX/2; 
Yvec = rMax-dX/2:-dX:-rMax+dX/2; 
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
[zMatrix] = zEval(rho,theta,maxDegree);
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
%
% Display map.
%
figure
imagesc(map,scale)
axis equal tight
axis off
colormap(palette)
colorbar
%
% Display min and max.
%
map = onMirror .* map;
minVal = min(min(map));
maxVal = max(max(map));
fprintf('\nfit min = %.3f\n',minVal);
fprintf('fit max = %.3f\n',maxVal);
