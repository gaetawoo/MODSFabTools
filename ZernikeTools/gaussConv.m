function smoothMap = gaussConv(map,mask,fwhmPix);
% smoothMap = convolve(map,mask,fwhmPix)
%
% Convolves map with Gaussian of width fwhmPix in pixels.
% Ignores pixels with mask = 0.
% Uses square kernel of width 2*fwhmPix.
%
% Define Gaussian kernel.
%
kernelCenter = round(fwhmPix)+1
kernelSize = 2*kernelCenter-1
gaussian = zeros(kernelSize);
xVec = -kernelCenter : kernelCenter; 
[x,y] = meshgrid(xVec,xVec);
argSq = (x.^2+y.^2)/fwhmPix^2;
gaussian = exp(-4*log(2)*argSq);
%
% Convolve map with Gaussian kernel.
%
[rows,cols] = size(map);
smoothMap = zeros(rows,cols);
for i=1:rows
    for j=1:cols
        if mask(i,j)
            sumWeightedMap = 0;
            sumWeight = 0;
            for iKer = 1 : kernelSize
                iSub = i+iKer-kernelCenter;
                if iSub>0 & iSub<=rows % Check that pixel is in array.
                    for jKer = 1 : kernelSize
                        jSub = j+jKer-kernelCenter;
                        if jSub>0 & jSub<=cols % Check that pixel is in array.
                            if mask(iSub,jSub)
                                weight = gaussian(iKer,jKer);
                                sumWeightedMap = sumWeightedMap...
                                    + weight*map(iSub,jSub);
                                sumWeight = sumWeight + weight;
                            end
                        end
                    end
                end
            end
            smoothMap(i,j) = sumWeightedMap/sumWeight;
        end
    end
end