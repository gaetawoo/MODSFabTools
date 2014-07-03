%%% Screen Clear and Variable Clear
clc;
clear;
 
%%% Linescan Adjustment (coeff = [linescan#,piston(+ is up),power(+ is convex)])
coeff = ...	
  [000,0.0,0.00; %% units in um
   005,0.4,0.00;
   010,0.4,0.00;
   015,0.5,0.00;
   020,-.1,0.00;
   025,-.2,0.00;
   030,0.0,0.00;
   035,0.3,0.07;
   040,-.4,-.04;
   045,-.5,0.02;
   050,-.7,0.00;
   055,-.7,0.00;
   060,-.5,0.00;
   065,-.7,0.00;
   070,-.8,0.00;
   075,0.3,0.05;
   080,0.4,0.00;
   085,0.4,0.04;
   090,0.6,0.02;
   095,0.0,-.04;
   100,0.0,0.00;
   105,-.6,0.00;
   110,-.8,0.00;
   115,0.6,0.00;
   120,-.6,0.00;
   125,0.0,0.00;
   130,0.0,0.00;
   135,0.0,0.00;
   140,-.7,0.00;
   145,-.7,0.00;
   150,-.9,0.00;
   155,-.8,0.00;
   160,-.9,0.00;
   165,-.7,0.00;
   170,-.9,0.00;
   175,0.2,0.05];
 
%%% Organizes Linescan Data
filename = '060511' %% filename is DATE_ANGLE.txt
 
k = 1;
B = [];
for i = 0:5:175
   %%% Load Linescan Data From File
   if i < 10 %% compensates numerical portion of filename for strokes less than 10
      filename2 = ['00' num2str(i)];	
   elseif i <100 & i > 5 %% compensates numerical portion of filename for less than 100
      filename2 = ['0' num2str(i)];
   else
      filename2 = num2str(i);
   end
 
   fid = fopen([filename '_' filename2 '.txt']);
   A = fscanf(fid,'%g   %g %g',[3 inf])';
   fclose(fid);
   
   A(:,1) = A(:,1)*1000; %% scale x-data from meters to mm
   A(:,2) = A(:,2)*1000; %% scale y-data from meters to mm
   A(:,3) = A(:,3)*1000000; %% scale z-data from meters to um
   
   %%% Adjust Linescan Data for Piston and Power
   for p = 1:size(A,1)
      A(p,3) = A(p,3)-(coeff(k,3)/10^5)*(7*(sqrt(A(p,1)^2 +A(p,2)^2))^2 - 1); %% power
   end
   A(:,3) = A(:,3) + coeff(k,2); %% piston
   
   %%% Write Out New Linescans After Piston and Power Adjustment
   Q = [];
   if A(1,1) == 0 & A(1,2) < 0 %% for x = 0 and y is negative
      for p = 1:size(A,1)
         if A(p,2) < 0
            Q(p,1) = -sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         else
            Q(p,1) = sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         end
      end
   elseif A(1,1) < 0 & A(1,2) < 0 %% for x is negative and y is negative
      for p = 1:size(A,1)	
         if A(p,2) < 0
            Q(p,1) = -sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         else
            Q(p,1) = sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         end
      end
   elseif A(1,1) < 0 & A(1,2) > 0 %% for x is negative and y is positive
      for p = 1:size(A,1)
         if A(p,2) > 0
            Q(p,1) = -sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         else
            Q(p,1) = sqrt(A(p,1).^2 + A(p,2).^2);
            Q(p,2) = A(p,3);
         end
      end
   end
   
   Q(:,1) = Q(:,1) - min(Q(:,1)); %% normalize
   dlmwrite([filename '_' filename2 '_rv3' '.txt'],Q,'-append','delimiter','\t'); %% write file
   
   %%% Plot New Linescan After Adjustment
   figure('Position',[6 65 1271 270])
   axes('DrawMode','fast','YGrid','on','FontWeight','bold') 
   axis([-500 500 -10 20])
   xlabel('Radial Position [mm]')
   ylabel('Surface Error [um]')
   title([filename '  ' filename2 ' Degrees'])
   box('on')
   hold('all')
   datasize = max(Q(:,1)) - min(Q(:,1));
   Q(:,1) = Q(:,1) - datasize/2;
   plot(Q(:,1),Q(:,2),'b')
 
   %%% Writes 3-Column Version of Linescan and 3-Column Version of All Linescans Combined
   A(:,1) = A(:,1)/1000; %% scale x-data from mm to meters	
   A(:,2) = A(:,2)/1000; %% scale y-data from mm to meters
   A(:,3) = A(:,3)/1000000; %% scale z-data from um to meters
   
   dlmwrite([filename '_' filename2 '_v3' '.txt'],A,'-append','delimiter','\t'); %% individual
   bigdatafile = ['MODSblue_' filename '_AllData_v3','.txt'];
   dlmwrite(bigdatafile,A,'-append','delimiter',','); %% all linescans
   
   B = [B;A];
   k = k + 1;
end
   
%%% Crop the Original Data or Not
cropornormal = 'normal'; %% 'crop' or 'normal' treatment of data
switch cropornormal
   case 'crop' %% makes a new matrix only with data within the specified annulus or aperture
      index = find(sqrt(B(:,1).^2 + B(:,2).^2) < 400 & sqrt(B(:,1).^2 + B(:,2).^2) > 76);
      C = zeros(size(index,1),3);
      for t = 1:size(index,1)
         C(t,1) = B(index(t),1);
         C(t,2) = B(index(t),2);
         C(t,3) = B(index(t),3);
      end
      x = C(:,1);
      y = C(:,2);
      z = C(:,3);
   case 'normal' %% does nothing to data
      x = B(:,1);
      y = B(:,2);
      z = B(:,3);
end
 
%%% Define Resolution for Interpolation (number of points in X and Y)
% n = 2*round((max(x)-min(x))/d); %% may require a lot of RAM and time
% n = 1000;  %% looks good
n = 400;  %% looks ok
% n = 20;  %% for fast processing only
 
%%% Make New X and Y Vectors Based on Resolution and Existing Data
X = linspace(min(x),max(x),n); %% interpolates new X dimention
Y = linspace(min(y),max(y),n); %% interpolates new Y dimention
scale = max([max(X) - min(X);max(Y) - min(Y)])/n; %% computes new pixel scale
 
%%% Interpolates Surface
M = griddata(x,y,z,X,Y','linear');
 
%%% Mask the New Array or Not
maskornot = 'not'; %% 'mask' or 'not' treatment of array
switch maskornot
   case 'mask' %% mask data array to specified clear aperture or annulus
      outer = 400; %% outer radial limit in mm
      inner = 76; %% inner radial limit in mm
      for i = 1:n
         for j = 1:n
            if sqrt((n/2 - i)^2 + (n/2 - j)^2) > outer/scale |...
                  sqrt((n/2 - i)^2 + (n/2 - j)^2) < inner/scale %% masks array data
               M(i,j) = NaN;
            end
         end
      end
   case 'not' %% do nothing to the array
end
 
%%% Create Border for Plot, Created Grid for Plot, Overlay Linescans on Surface Map
border = 'yes'; %% 'yes' or 'no' to create border for plot
grid = 'no'; %% 'yes' or 'no' to create grid for plot
linescanoverlay = 'no'; %% 'yes' or 'no' to overlay linescan data on surface map
 
switch border
   case 'yes' %% creates border of specified annulus or clear aperture
      outer = 400; %% units in mm	
      inner = 76; %% units in mm
      y1 = -outer:.1:outer;
      y2 = -inner:.1:inner;
      x1 = [sqrt(outer^2 - y1.^2),-sqrt(outer^2 - y1.^2)];
      x2 = [sqrt(inner^2 - y2.^2),-sqrt(inner^2 - y2.^2)];
   case 'no' %% makes no border
end
 
switch grid
   case 'yes'
      gridd = -420:420; %% units in mm
      gridd2 = -420:70:420; %% units in mm
   case 'no'
end
 
%%% Normalize and Correct for Any NaN Values
M = M - min(min(M)) + eps; %% normalize data
M(isnan(M)==1) = 0; %% correct NaN values so computations and plotting don't crash
 
%%% Set Colormap and Smooth Surface (removes interpolation artifacts)
colormapPV = [min(min(M))+0,min(min(M))+11]; %% set colormap to defined scale
M = Array_Smooth(M,830,0,415,10,'jet',colormapPV,colormapPV,X,Y); %% B.Martin program
                                                                  %% 10 mm FWHM smoothing
%%% Plot Surface
figure('Position',[350 100 900 800])
pcolor(X,Y,M) %% plots array
shading interp
colormap jet
caxis(colormapPV) %% define colorscale
colorbar(...
  'Box','on',...
  'YMinorTick','on',...
  'YTick',min(colormapPV):1:max(colormapPV),...
  'YLim',[colormapPV]);
title(['MODS Blue CMM Data Surface Map - ' filename])
axis equal
axis tight
xlabel('X [mm]')
ylabel('Y [mm]')
hold all
 
switch border %% plots border
   case 'yes'
      plot(x1,[y1,y1],'k','LineStyle','none','Marker','.','MarkerSize',1)
      plot(x2,[y2,y2],'k','LineStyle','none','Marker','.','MarkerSize',1)
   case 'no'
end
 
switch grid %% plots grid
   case 'yes'
      for w = 1:size(gridd2,2)
         plot(gridd,gridd2(w),'k')
         plot(gridd2(w),gridd,'k')
      end
   case 'no'
end
 
switch linescanoverlay %% plots linescans over surface data to show were real data exists
   case 'yes'
      plot3(x,y,z,'r.','MarkerSize',2)
   case 'no'
end
 
%%% Write Out Surface Map in n by n Array
dlmwrite(['MODSblue_' filename '_array830.txt'],M,'-append','delimiter',' '); %% save file
 
%%% Write Out Surface Map in 3-Column Matrix
data_array = zeros(n^2,3); %% populate matrix for faster processing
 
k = 1;
for i = 1:n
   for j = 1:n
      data_array(k,1) = X(j); %% transform array to x-coordinate
      data_array(k,2) = Y(i); %% transform array to y-coordinate
      data_array(k,3) = M((i),(j)); %% transform array to z-coordinate
      k = k + 1;
   end
end
 
indexx = find(data_array(:,3) ~= 0 |...
   sqrt(data_array(:,1).^2 + data_array(:,2).^2) < 400); %% only save points where surface is
 
ohhs = zeros(size(indexx,1),3); %% populate matrix for faster processing
for i = 1:size(indexx,1)
   ohhs(i,1) = data_array(indexx(i),1); %% populate column with surface data in x
   ohhs(i,2) = data_array(indexx(i),2); %% populate column with surface data in y
   ohhs(i,3) = data_array(indexx(i),3); %% populate column with surface data in z
end
 dlmwrite(['MODSblue_' filename '_xyz100.txt'],ohhs,'-append','delimiter',' '); %% save file
 
%%% Plot Surface Data as 3-D Point Map
cold_fusion(ohhs,colormap) %% B.Steward program to visualize 3-D points with custom colormap
