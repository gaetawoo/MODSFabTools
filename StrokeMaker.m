%%% Screen Clear and Variable Clear
clc;
clear;
 
%%% Load Data (file format is 2-column matrix, space delimited)
filename1 = 'MODSred_';	
filename2 = '071025';
filename3 = '_ARP';
filename = [filename1 filename2 filename3]; %% filename is PROJECT_DATE_DATATYPE.txt
fileext = '.txt';
[F1x,F1y] = textread([filename fileext],'%f %f','delimiter',' '); %% units of [mm,um]
totalsize = size(F1x,1);
 
%%% Make a Changable Version and an Original Version (F1y2 is original version)
F1y2 = F1y;
 
%%% Stroke Parameters
lapsize = 4*25.4; %% units of mm
partradius = 420; %% units of mm
 
%%% User-Input Duration of Each Stroke
duration = ... %% duration of each stroke in hours [strokenumber,hours]
  [1,3;
   2,1;
   3,1;
   4,1;
   5,2;
   6,1;
   7,3;
   8,1;
   9,3;
  10,1;
  11,4;
  12,1;
  13,4;
  14,1;
  15,4;
  16,1;
  17,3;
  18,1;
  19,3;
  20,1;
  21,3;
  22,1;
  23,1;
  24,1;
  25,1;
  26,1;
  27,8]*5; %% scaling factor for quickly scaling entire run
 
%%% Plot Parameters
figure('Position',[6 65 1271 270],'NumberTitle','off','Name',filename)
axes1 = axes(...
  'DrawMode','fast',...
  'FontWeight','bold',...
  'YGrid','on');
axis(axes1,[0 420 -5 15])
xlabel('Radial Position [mm]')
ylabel('Surface Error or Removal [um]')
title(['MODS Red ' datestr(now,'yymmdd') ' Scanning Run'])
box('on')
hold('all')
 
%%% Prepare Vector and Matrix (pre-populates for faster processing)
removal_center = zeros(size(duration,1),1);
bins = zeros(size(duration,1),6);
for p = 1:size(duration,1)
   bins(p,1) = p;
   bins(p,2) = duration(p,2);
end
 
%%% Processing Loop (surface profile is adjusted by stroke removal one stroke at a time)
for i = 1:size(duration,1)
   
   %%% Load Stroke Removal Data
   filename4 = '71019_'; %% removal profile filename is DATE_## (2 digit stroke number)
   if i < 10 %% compensates numerical portion of filename for strokes less than 10
      filenum = ['0' num2str(i)];
   else
      filenum = num2str(i);
   end
   
   %%% Adjust Surface Profile by Stroke Removal
   [x_rev,y_rev] = textread([filename4 filenum],'%f %f','delimiter',' '); %% units of [%,um]
   index = max(find(y_rev == min(y_rev))); %% find data index of maximum removal for stroke
   center_rev = x_rev(index); %% get location on part of maximum removal for stroke
   x_rev = x_rev*partradius; %% convert x-axis from % of part to mm	
   sep = (0:2:F1x(size(F1x,1)))'; %% make x-axis data points for new interpolated removal
   y_rev_new = spline(x_rev,y_rev,sep); %% new interpolated removal in um
   F1y = F1y + y_rev_new*duration(i,2); %% adjust surface profile based on stroke removal
   indv_lap = y_rev_new*duration(i,2); %% individual stroke removal
   removal_center(i) = center_rev*partradius; %% location of max stroke removal in mm
   notch = min(indv_lap):.01:0; %% magnitude of maximum removal (negative)
   
   %%% Loading Stroke Parameter Data
   [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19]...
   = textread([filename4 filenum '.TXT'],...
   '%s %s %f %f %f %f %f %f %f %f %f %f %f %s %s %f %f %f %f','delimiter',',');
   bins(i,3) = A3*partradius - lapsize/2; %% converts inner lap pin location to inner lap edge
   bins(i,4) = A4*partradius + lapsize/2; %% converts outer pin location to outer edge
   bins(i,5) = ((A4 + A3)/2)*partradius; %% calculates center of stroke location
   bins(i,6) = bins(i,4) - bins(i,3) - lapsize; %% calculates stroke size
      
   plot(F1x,indv_lap,'m',removal_center(i),notch,'k.',bins(i,5),notch,'w');
   hold on
end
 
%%% Final Plot
plot(F1x,F1y2,'b') %% plots original profile
plot(F1x,F1y,'g','LineWidth',2) %% plots new profile after removal
plot(F1x,(F1y-F1y2),'r') %% plots removal profile
plot(75:.1:397,(max(F1y(37:200))),'r') %% plots line showing new max surface height
plot(75:.1:397,(min(F1y(37:200))),'r') %% plots line showing new min surface height
plot(397,min(F1y(37:200)):.01:max(F1y(37:200)),'r') %% plots line showing outer clear aperture
plot(75,min(F1y(37:200)):.01:max(F1y(37:200)),'r') %% plots line showing inner clear aperture
hold off
 	
%%% Write Calculated Stroke Parameters to Excel
xlswrite([filename '_ScanningRun.xls'],bins,'Scanning Parameters')
 
%%% Display Stroke Parameters on Screen
for q = 1:size(bins,1) %% converts stroke parameters into inches and rounds to tenths place
bins(q,3) = (round(10*(bins(q,3)/25.4)))/10;
bins(q,4) = (round(10*(bins(q,4)/25.4)))/10;
bins(q,5) = (round(10*(bins(q,5)/25.4)))/10;
bins(q,6) = (round(10*(bins(q,6)/25.4)))/10;
end
 
'       #    |Hours | Inner Edge|Outer Edge| Center  | Stroke Size' %% Column Headings
bins %% screen print of stroke parameters
totalTIME = sum(bins(:,2)) %% total run time calculated by summing durations of each stroke
 	
stats(1,1) = 1; %% row heading for original profile
stats(2,1) = 2; %% row heading for new profile
stats(1,2) = max(F1y2(37:200))-min(F1y2(37:200)); %% peak-to-valley of original profile
stats(1,3) = sqrt(mean(F1y2(37:200).^2)-(mean(F1y2(37:200)))^2); %% rms of original profile
stats(2,2) = max(F1y(37:200))-min(F1y(37:200)); %% peak-to-valley of new profile
stats(2,3) = sqrt(mean(F1y(37:200).^2)-(mean(F1y(37:200)))^2); %% rms of new profile
 
stats
