function ionFigure(file,D,rMin,rMax,fwhm,palette,beforeScale,afterScale,maskVal,boundaryVal,delta)
%
% ionFigure(file,D,rMin,rMax,fwhm,palette,beforeScale,afterScalemaskVal,boundaryVal,delta)
%
% Performs synthetic ion figuring of NST mirror
%
% INPUT
% file      name of .int file, with full path and extension
% D         diameter of array in .int file
% rMin      of subaperture to be processed and displayed
% rMax      of subaperture to be processed and displayed
% fwhm      of Gaussian ion beam
% palette   Matlab colormap (string)
% beforeScale	[min,max] for color map for raw and smoothed data
% afterScale	[min,max] for color map for difference (residual data)
% maskVal   value to display for region outside subaperture
% boundaryVal   value to display for boundary of subaperture
% delta     half-width of boundary

%
% Read map. loadCodeV() returns wavefront, so divide by 2 for surface. 
%
[map,mask,rows,cols] = loadCodeV(file,0);
map = map / 2;
%
% Use dimensions of unit circle to assign x and y coordinates to pixels.
%
UCdiam = rows;
dX = D/UCdiam;
fwhmPix = fwhm/dX;
Xvec = -D/2+dX/2 : dX : D/2-dX/2; 
if length(Xvec) ~= cols
    error('Error assigning coordinates to pixels.')
end
[X,Y] = meshgrid(Xvec,Xvec);
R = sqrt(X.*X+Y.*Y);
onMirror = (R>=rMin & R<=rMax);
ptsOnMirror = sum(sum(onMirror))
map = onMirror .* map;
mask = onMirror .* mask;
%
% Define boundary of clear aperture.
%
innerICA = rMin-delta;
outerICA = rMin+delta;
innerOCA = rMax-delta;
outerOCA = rMax+delta;
boundary = (R>innerICA & R<outerICA) | (R>innerOCA & R<outerOCA);
%
% Shrink array to unmasked area.
%
% imin = rows;
% imax = 1;
% jmin = cols;
% jmax = 1;
% for i=1:rows
%     for j=1:cols
%         if mask(i,j)
%             if i<imin
%                 imin = i;
%             elseif i>imax
%                     imax = i;
%             end
%             if j<jmin
%                 jmin = j;
%             elseif j>jmax
%                 jmax = j;
%             end
%         end
%     end
% end
% map = map(imin:imax,jmin:jmax);
% [rows,cols] = size(map)
% mask = mask(imin:imax,jmin:jmax);
ptsInMask = sum(sum(mask))
mean = sum(sum(map))/ptsInMask
map = map - mean;
map = mask .* map;
rms = sqrt(sum(sum(map.*map))/ptsInMask);
maxVal = max(max(map));
minVal = min(min(map));
%
% Convolve map with Gaussian.
%
smoothMap = gaussConv(map,mask,fwhmPix);
meanSmooth = sum(sum(smoothMap))/ptsInMask
smoothMap = smoothMap - meanSmooth;
smoothMap = mask .* smoothMap;
rmsSmooth = sqrt(sum(sum(smoothMap.*smoothMap))/ptsInMask);
maxSmooth = max(max(smoothMap));
minSmooth = min(min(smoothMap));
%
% Subtract smoothed map from original map.
%
residMap = map-smoothMap;
meanResid = sum(sum(residMap))/ptsInMask
% residMap = residMap - meanResid;
rmsResid = sqrt(sum(sum(residMap.*residMap))/ptsInMask);
maxResid = max(max(residMap));
minResid = min(min(residMap));
%
% Assign values to masked areas and boundary.
%
for i=1:rows
    for j=1:cols
        if ~mask(i,j)
            map(i,j) = maskVal;
            smoothMap(i,j) = maskVal;
            residMap(i,j) = maskVal;
        end
        if boundary(i,j)
            map(i,j) = boundaryVal;
            smoothMap(i,j) = boundaryVal;
            residMap(i,j) = boundaryVal;
        end
    end
end
%
% Display maps.
%
close all
figure
imagesc(map,beforeScale)
axis equal tight
axis off
colormap(palette)
colorbar
title(['original: ',num2str(rms,'%.1f'),' nm rms, range = (',...
    num2str(minVal,'%.0f'),',',num2str(maxVal,'%.0f'),')'])
figure
imagesc(smoothMap,beforeScale)
axis equal tight
axis off
colormap(palette)
colorbar
title(['smoothed: ',num2str(rmsSmooth,'%.1f'),' nm rms, range = (',...
    num2str(minSmooth,'%.0f'),',',num2str(maxSmooth,'%.0f'),')'])
figure
imagesc(residMap,afterScale)
axis equal tight
axis off
colormap(palette)
colorbar
title(['residual: ',num2str(rmsResid,'%.1f'),' nm rms, range = (',...
    num2str(minResid,'%.0f'),',',num2str(maxResid,'%.0f'),')'])