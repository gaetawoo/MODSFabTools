%function subapLOTIS(file,lclp,uclp)


for sss = 1:1
clear D angles data_array i map newN rows windy zCoef Dhole datareport index mask q scale x zMatrix RMS0400 aptnum diam j nVec r t0 x_subap RMS1000 cols elVec k newData r_subap uclp y_subap RMS2500 compression file lclp newData2 rho_a windx z subData subData2 windy2 windx2;
datestr(now)
t0 = clock;
lclp = -75;
uclp = +75;
compression = 1;
tic
switch sss
case 1
filename = '060121'
case 2
filename = '060122'
case 3
filename = '060123'
case 4
filename = '051118'
case 5
filename = '051115'
case 6
filename = '051109'
case 7
filename = '051107'
case 8
filename = '051031'
case 9
filename = '051019'
case 10
filename = '051018'
case 11
filename = '051014'
case 12
filename = '051011'
case 13
filename = '051007'
case 14
filename = '051005'
case 15
filename = '051003'
end
file = ['\primary data\' filename '.int'];
%
% displayInt(file,rmax,rmin,lclp,uclp)
%
% Displays CodeV file at full resolution.
%
% INPUT
% file      file name with full path and extension
% rmax      mask out area outside this radius (m)
% rmin      mask out area inside this radius (m)
% lclp      value corresponding to bottom of colorbar
% uclp      value corresponding to top of colorbar
%
D = 6.46; % meters
Dhole = 0.942; % meters
%
% Read map. loadCodeV() returns wavefront, so divide by 2 for surface. 
% In order to look right after reading with loadCodeV(), LBT maps need to 
% be flipped around y axis.
%
[map,mask,rows,cols] = loadCodeV(file,0);
% map = fliplr(map);
map = map / 2;

% Display map.
%
% close all
% figure(4)
% imagesc(map,[lclp uclp])
% axis equal tight
% axis off
% colormap('jet')
% colorbar

total_rows = (rows/compression)*(cols/compression);
data_array = zeros(total_rows,3);

k = 1;
for i = 1:compression:rows
	for j = 1:compression:cols
		data_array(k,1) = j;
		data_array(k,2) = 441-i;
		data_array(k,3) = map((i),(j));
		k = k + 1;
	end
end
% data_array
i = 1;
j = 1;
k = 1;
toc

% 6.5 m aperture-----------------------------------------------------------
x_subap = (440/compression * 0.5)*compression;
y_subap = (440/compression * 0.5)*compression;
r_subap = (440/compression)/2;
r = sqrt((data_array(:,1)-x_subap).^2 + (data_array(:,2)-y_subap).^2);


scale = 6500/440;
index = find(r <= r_subap*.9938 & r >= r_subap*.1449);
newN = length(index);
newData = zeros(newN,3);
for i = 1:newN
	newData(i,1) = data_array(index(i),1);
	newData(i,2) = data_array(index(i),2);
	newData(i,3) = data_array(index(i),3);
end

normaliz_rho = (((6500*.9938)/2)/scale);
rho = sqrt((newData(:,1)-x_subap).^2 + (newData(:,2)-y_subap).^2)/normaliz_rho;
theta = atan2(((newData(:,2)-y_subap)),((newData(:,1)-x_subap)));
z = newData(:,3);

[zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,z,1);
newData2(:,1) = newData(:,1);
newData2(:,2) = newData(:,2);
newData2(:,3) = newData(:,3) - zMatrix*zCoef;
% newData

datareport(1,1) = sqrt(mean(newData2(:,3).^2)-(mean(newData2(:,3)))^2);
datareport(2,1) = sqrt(mean(newData2(:,3).^2)-(mean(newData2(:,3)))^2);
datareport(3,1) = sqrt(mean(newData2(:,3).^2)-(mean(newData2(:,3)))^2);

RMS6500(k,3) = sqrt(mean(newData2(:,3).^2)-(mean(newData2(:,3)))^2);
RMS6500(k,1) = 0;
RMS6500(k,2) = 0;
RMS6500(k,4) = max(newData2(:,3)) - min(newData2(:,3));

% 2.5 m aperture-----------------------------------------------------------

rho_a = [1720.925/scale;1979.85/scale];  % select radial zones to sample from
diam = 2500;	% diameter of supaperture
normaliz_rho = ((diam/2)/scale);

for i = 1:size(rho_a,1)
	aptnum = rho_a(i)*2*pi()/(diam/scale);	% determines number of subapertures that can fit in the radial position
	switch mod(ceil(aptnum),2)	% rounds # of subapertures
		case 1
		aptnum = ceil(aptnum) + 1;
		case 0
		aptnum = ceil(aptnum);
	end
	angles = 360/aptnum;	% detemines the angle seperation of each aperture within a radial position
	x = (0:1:aptnum-1)';	% used only here to induce a 30 degree offset from where to start the first subap.
	switch i
		case 1
			angles = angles*x+30;
		case 2
			angles = angles*x;
	end	
	for j = 1:1:size(angles,1)	
		windx = rho_a(i)*cos(angles(j)*pi()/180)+x_subap;   % gives the position of the center of the subap
		windy = rho_a(i)*sin(angles(j)*pi()/180)+y_subap;
		windx2 = rho_a(i)*cos(angles(j)*pi()/180)*scale;	% gives physical x,y coordinates of center of subap
		windy2 = rho_a(i)*sin(angles(j)*pi()/180)*scale;
	
		r = sqrt((newData(:,1)-windx).^2 + (newData(:,2)-windy).^2);	% creates vector of radial positions
		index = find(r <= ((diam/2)/scale));	% finds the indices of the radial points within the window of the subap
		newN = length(index);
		subData = zeros(newN,3);
		for q = 1:newN	% creates new matrix using only data within the window
			subData(q,1) = newData(index(q),1);
			subData(q,2) = newData(index(q),2);
			subData(q,3) = newData(index(q),3);
		end
		rho = sqrt((subData(:,1)-windx).^2 + (subData(:,2)-windy).^2)/normaliz_rho; 
		theta = atan2(((subData(:,2)-windy)),((subData(:,1)-windx)));
		z = subData(:,3);

		[zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,z,1); % send out rho,theta,z and receive back detiled zernike fit for subap
		subData2(:,1) = subData(:,1);
		subData2(:,2) = subData(:,2);
		subData2(:,3) = subData(:,3) - zMatrix*zCoef;
		% newData

		RMS2500(k,1) = windx2;
		RMS2500(k,2) = windy2;
		RMS2500(k,3) = sqrt(mean(subData2(:,3).^2)-(mean(subData2(:,3)))^2);		
		RMS2500(k,4) = max(subData2(:,3)) - min(subData2(:,3));
		
		k = k + 1;
		clear subData subData2 rho theta
	end
end

k-1;
datareport(2,2) = max(RMS2500(:,3));
datareport(3,2) = min(RMS2500(:,3));
datareport(1,2) = mean(RMS2500(:,3));

clear rho_a diam
i = 1;
j = 1;
k = 1;

% 1.0 m aperture-----------------------------------------------------------

rho_a = [970.925/scale;1410.65625/scale;1850.3875/scale;2290.11875/scale;2729.85/scale];
diam = 1000;
normaliz_rho = ((diam/2)/scale);

for i = 1:size(rho_a,1)
	aptnum = rho_a(i)*2*pi()/(diam/scale);
	switch mod(ceil(aptnum),2)
		case 1
		aptnum = ceil(aptnum) + 1;
		case 0
		aptnum = ceil(aptnum);
	end
	angles = 360/aptnum;
	x = (0:1:aptnum-1)';
	angles = angles*x;
	for j = 1:1:size(angles,1)
		windx = rho_a(i)*cos(angles(j)*pi()/180)+x_subap;
		windy = rho_a(i)*sin(angles(j)*pi()/180)+y_subap;
		windx2 = rho_a(i)*cos(angles(j)*pi()/180)*scale;
		windy2 = rho_a(i)*sin(angles(j)*pi()/180)*scale;
	
		r = sqrt((newData(:,1)-windx).^2 + (newData(:,2)-windy).^2);
		index = find(r <= ((diam/2)/scale));
		newN = length(index);
		subData = zeros(newN,3);
		for q = 1:newN
			subData(q,1) = newData(index(q),1);
			subData(q,2) = newData(index(q),2);
			subData(q,3) = newData(index(q),3);
		end
		rho = sqrt((subData(:,1)-windx).^2 + (subData(:,2)-windy).^2)/normaliz_rho;
		theta = atan2(((subData(:,2)-windy)),((subData(:,1)-windx)));
		z = subData(:,3);

		[zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,z,1);
		subData2(:,1) = subData(:,1);
		subData2(:,2) = subData(:,2);
		subData2(:,3) = subData(:,3) - zMatrix*zCoef;
		% newData

		RMS1000(k,3) = sqrt(mean(subData2(:,3).^2)-(mean(subData2(:,3)))^2);
		RMS1000(k,1) = windx2;
		RMS1000(k,2) = windy2;		
		RMS1000(k,4) = max(subData2(:,3)) - min(subData2(:,3));
		k = k + 1;
		clear subData subData2 rho theta
	end
end

k-1;
datareport(2,3) = max(RMS1000(:,3));
datareport(3,3) = min(RMS1000(:,3));
datareport(1,3) = mean(RMS1000(:,3));

clear rho_a diam
i = 1;
j = 1;
k = 1;

% 0.4 m aperture-----------------------------------------------------------

rho_a = [670.925/scale;933.02778/scale;1195.1305556/scale;1457.23333/scale;1719.33611/scale;1981.438889/scale;2243.541667/scale;2505.6444/scale;2767.7472222/scale;3029.85/scale];
diam = 400;
normaliz_rho = ((diam/2)/scale);

for i = 1:size(rho_a,1)
	aptnum = rho_a(i)*2*pi()/(diam/scale);
	switch mod(ceil(aptnum),2)
		case 1
		aptnum = ceil(aptnum) + 1;
		case 0
		aptnum = ceil(aptnum);
	end
	angles = 360/aptnum;
	x = (0:1:aptnum-1)';
	angles = angles*x;
	for j = 1:1:size(angles,1)
		windx = rho_a(i)*cos(angles(j)*pi()/180)+x_subap;
		windy = rho_a(i)*sin(angles(j)*pi()/180)+y_subap;
		windx2 = rho_a(i)*cos(angles(j)*pi()/180)*scale;
		windy2 = rho_a(i)*sin(angles(j)*pi()/180)*scale;

		
		r = sqrt((newData(:,1)-windx).^2 + (newData(:,2)-windy).^2);
		index = find(r <= ((diam/2)/scale));
		newN = length(index);
		subData = zeros(newN,3);
		for q = 1:newN
			subData(q,1) = newData(index(q),1);
			subData(q,2) = newData(index(q),2);
			subData(q,3) = newData(index(q),3);
		end
		rho = sqrt((subData(:,1)-windx).^2 + (subData(:,2)-windy).^2)/normaliz_rho;
		theta = atan2(((subData(:,2)-windy)),((subData(:,1)-windx)));
		z = subData(:,3);

		[zCoef,zMatrix,nVec,elVec] = zFit(rho,theta,z,1);
		subData2(:,1) = subData(:,1);
		subData2(:,2) = subData(:,2);
		subData2(:,3) = subData(:,3) - zMatrix*zCoef;
		% newData

		RMS0400(k,3) = sqrt(mean(subData2(:,3).^2)-(mean(subData2(:,3)))^2);
		RMS0400(k,1) = windx2;
		RMS0400(k,2) = windy2;
		RMS0400(k,4) = max(subData2(:,3)) - min(subData2(:,3));
		k = k + 1;
		clear subData subData2 rho theta
	end
end

k-1;
datareport(2,4) = max(RMS0400(:,3));
datareport(3,4) = min(RMS0400(:,3));
datareport(1,4) = mean(RMS0400(:,3));
%--------------------------------------------------------------------------
datareport
xlswrite([filename '_SubAp_Analysis.xls'],datareport,'Data Report')
xlswrite([filename '_SubAp_Analysis.xls'],RMS6500,'6.5m Aperture')
xlswrite([filename '_SubAp_Analysis.xls'],RMS2500,'2.5m Aperture')
xlswrite([filename '_SubAp_Analysis.xls'],RMS1000,'1.0m Aperture')
xlswrite([filename '_SubAp_Analysis.xls'],RMS0400,'0.4m Aperture')
fprintf('This calculation took %g seconds to run.\n',(etime(clock,t0)))
end
