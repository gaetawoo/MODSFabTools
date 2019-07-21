function [wavefront,mask,rows,cols]=loadCodeV(file_name,with_graph)
%function [wavefront,mask,rows,cols]=loadCodeV(file_name,with_graph)
%----------------------------------------------------------------------------
% this function accepts the name of a file containing an image in codeV format
% and returns this image as a 2D array
%
% INPUTS
% file_name				the file name of the codeV coded image
% with_graph			flag to enable or disable the graphing of the results
%
% OUTPUTS
% wavefront				a rows x cols 2D array containing the wavefront in nm
% mask					pupil mask: 
%							a 2D table rows  x cols in size, which has a value of 1 
% 							inside pupil and 0 outside.
%
% HISTORY
% V0.1 Created: 	by F. Wildi, March 99, Steward Observatory
% V0.2 Modified: 	by F. Wildi, June 99, to support running with 'n_sample' input
% V0.3 Modified, 29.7.99: 	by F. Wildi to be more robust to different sizes of images and 
%						to remove the transpose perticular to the original WFS calibration 
%						problem
% V0.4 Modified 14.9.01 by B. Martin, to remove resampling and return wavefront rather
%						than surface. Name changed from codeV2Matlab to loadCodeV
%----------------------------------------------------------------------------

%
% Read the file
%
filein=fopen(file_name);
line1=fread(filein,'char');
length_data=length(line1);
fclose(filein);

ind_num=0;  % numbers the numerical values in the end array
ind_line=1; % numbers the characters in the input array
end_while=0;
start=1;
%
% Extract the header, reformatting it and saving it to a temporary file 
%
comment_array=[];
while (ind_line<length_data) && ~end_while
   ind_line=ind_line+1;
   % extraction
   if sum(line1(ind_line-1:ind_line)==[13;10])==2  && ...
         ~strcmp(char(line1(ind_line+1:ind_line+3))','GRD') % eol in the comment header
      ind_num=ind_num+1;
      comment_array=char(comment_array,line1(start:ind_line)');
      start=ind_line+1;
   elseif strcmp(char(line1(ind_line+1:ind_line+3))','GRD') % beginning factors line
      ind_num=ind_num+1;
      comment_array=char(comment_array,line1(start:ind_line)');
      comment_array=[char(ones(size(comment_array,1),1)*'%') comment_array];
      ind_line=ind_line+1;
      start=ind_line;
      while (ind_line<length_data) && ~(sum(line1(ind_line-1:ind_line)==[13;10])==2)
         ind_line=ind_line+1;
      end
      comment_array=char(comment_array,line1(start:ind_line)');
      start=ind_line+1;
      % reformatting/saving
      comment_file='temp_header.txt';
      fid=fopen(comment_file,'w');
      n_lines=size(comment_array,1);
      fwrite(fid,comment_array(n_lines,:),'uchar');
      fprintf(fid,'\n');
      for ind=1:n_lines-1
         %save(comment_file,'comment_array','-ascii');
         fwrite(fid,comment_array(ind,:),'uchar');
         fprintf(fid,'\n');
      end
      fclose(fid);
      %save(comment_file,'comment_page','-ascii')
      data_file='temp_data.txt';
      data=line1(start:length(line1));
      %save(data_file,'data','-ascii')
      fid=fopen(data_file,'w');
      fprintf(fid,'%s',data);
      fclose(fid);
      end_while=1;
   end
end
%
% Loading the temporary file containing the header
%
fid=fopen(comment_file,'r');
%
% Extact certain header data
%
coding_type			=fscanf(fid,'%3s',1);
rows				=fscanf(fid,'%d',1);
cols				=fscanf(fid,'%d',1);
data_type			=fscanf(fid,'%3s',1);
wavelength_symbol	=fscanf(fid,'%3s',1);
wavelength			=fscanf(fid,'%f',1);
scale_symbol		=fscanf(fid,'%3s',1);
scale				=fscanf(fid,'%d',1);
invalid_symbol		=fscanf(fid,'%3s',1);
invalid_value		=fscanf(fid,'%d',1);
fclose(fid);
qq=sprintf('coding:%3s rows:%3d cols:%3d data_type:%3s lambda:%f %3s:%3d %3s:%3d',...
   coding_type,rows,cols,data_type,wavelength,scale_symbol,scale,invalid_symbol,...
   invalid_value);
disp(qq)
%
% Load the temporary file containing the data
%
fid=fopen(data_file,'r');
raw_data=fscanf(fid,'%d',rows*cols);
%
% reshape() below fills up matrix data column by column. We want it done row by row.
%
data = reshape(raw_data,rows,cols);
data = data';
%
% Display the raw image
%
%close all
if with_graph
   figure
   imagesc(data)
   axis equal tight
   colorbar
   xlabel('x axis')
   ylabel('y axis')
   qq=sprintf('raw image');
   title(qq)
end
%
% Define pupil
%
mask = (data ~= invalid_value);
%
% Scale the data to waves
%
data=data/scale;
%
% Convert to wavefront if not already there.
%
if ~strcmp(data_type,'WFR')
   data=data*2;
end
%
% Convert to nm (the wavelength is in micron)
%
data=data*wavelength*1000;
%
% Apply pupil mask, subtract mean, and reapply mask
%
wavefront = mask .* data;
validData = sum(sum(mask));
wavefront = wavefront - sum(sum(wavefront))/validData;
wavefront = mask .* wavefront;
if with_graph
   figure
   imagesc(wavefront);
   axis equal tight
   colorbar
   xlabel('x axis')
   ylabel('y axis')
   qq=sprintf('image with pupil mask');
   title(qq)
end
