%%% Program and GUI Initialization
function varargout = LinescanAlignment(varargin)
gui_Singleton = 1; %% Begin initialization code - DO NOT EDIT
gui_State = struct('gui_Name',       mfilename, ...
   'gui_Singleton',  gui_Singleton, ...
   'gui_OpeningFcn', @LinescanAlignment_OpeningFcn, ...
   'gui_OutputFcn',  @LinescanAlignment_OutputFcn, ...
   'gui_LayoutFcn',  [] , ...
   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
   [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end %% End initialization code - DO NOT EDIT
 
%%% Executes Just Before GUI is Made Visible
function LinescanAlignment_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject; %% Choose default command line output for LinescanAlignment
guidata(hObject, handles); %% Update handles structure
if nargin == 3,
   initial_dir = pwd;
elseif nargin > 4
   if strcmpi(varargin{1},'dir')
      if exist(varargin{2},'dir')
         initial_dir = varargin{2};
      else
         errordlg('Input argument must be a valid directory','Input Argument Error!')
         return
      end
   else
      errordlg('Unrecognized input argument','Input Argument Error!');
      return;
   end
end
load_listbox(initial_dir,handles) %% Populate the listbox
 
%%% Outputs of this Function are Returned to the Command Line
function varargout = LinescanAlignment_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
 
%%% Loads Linescan Data and Logfile Data On a Double Click in the Listbox
function varargout = listbox1_Callback(hObject, eventdata, handles)
get(handles.figure1,'SelectionType');
if strcmp(get(handles.figure1,'SelectionType'),'open')
   index_selected = get(handles.listbox1,'Value');
   file_list = get(handles.listbox1,'String');
   filename = file_list{index_selected};
   if  handles.is_dir(handles.sorted_index(index_selected))
      cd (filename)
      load_listbox('C:\MODS',handles)
   else
      if get(handles.radiobutton1,'Value') == 0 & get(handles.radiobutton2,'Value') == 0
         errordlg('Please choose a scan method.','No Scan Method Selected','modal')
      else
         [path,name,ext,ver] = fileparts(filename);
         set(handles.togglebutton2,'Value',0,'String','Plot with Cut/Clear Aperture Lines');
         switch ext
            case '.DAT'
               handles.file1 = name;
               file2 = [name,ext];
               
               set(handles.text22,'String',file2);
               
               [F1x, F1y, F1z, F1i, F1j, F1k] = textread(file2,... %% loads file into matrix
                  '%f %f %f %f %f %f','delimiter',',');
               totalsize = size(F1x,1);
               
               handles.vectsize = totalsize - 6;
               X = F1x(7:totalsize);
               Y = F1y(7:totalsize);
               Z = F1z(7:totalsize);
               Xv = F1i(7:totalsize);
               Yv = F1j(7:totalsize);
               Zv = F1k(7:totalsize);
               handles.probedia = F1x(1);
               
               if get(handles.radiobutton1,'Value')
                  for iii = 1:size(X,1) %% compensate for probe
                     handles.X(iii,1) = X(iii) - Xv(iii)*(handles.probedia/2);
                     handles.Y(iii,1) = Y(iii) - Yv(iii)*(handles.probedia/2);
                     ZZ(iii,1) = Z(iii) - Zv(iii)*(handles.probedia/2);
                  end
               else
                  handles.X = X;
                  handles.Y = Y;
                  ZZ = Z;
               end
               
               set(handles.text19,'String',num2str(handles.probedia));
               
               handles.piston = 0;
               handles.Z = ZZ - max(ZZ);
               set(handles.togglebutton3,'Value',0);
               set(handles.togglebutton3,'String','Plot OPD Surface Data')
               set(handles.togglebutton4,'Value',0);
               set(handles.togglebutton4,'String','Plot Earlier Line')
               handles.hflip = 0;
               set(handles.togglebutton1,'Value',0);
               
               numeros = size(handles.file1,2);
               handles.angle = str2num(handles.file1(numeros-2:numeros));
                              
               set(handles.edit1,'String',handles.file1(numeros-2:numeros));
               
               try
                  fid_a = fopen('linescanlog.txt','r'); %% load and read logfile
                  Q = textscan(fid_a,'%50c','delimiter',',');
                  s = Q{:};
                  n_out = 0;
                  for n = 1:size(s,1)
                     if findstr(name,s(n,:))
                        n_out = n;
                        dataline = s(n_out,:);
                        [o,o,XC,YC,hflip,angle,saved] = strread(dataline,...
                           '%s %s %s %s %s %s %s','delimiter',',');
                        
                        handles.XC = str2num(XC{:});
                        handles.YC = str2num(YC{:});
                        handles.hflip = str2num(hflip{:});
                        handles.angle = str2num(angle{:});
                        handles.saved = str2num(saved{:});
                        
                        set(handles.togglebutton1,'Value',handles.hflip);
                        set(handles.edit1,'String',num2str(handles.angle));
                        set(handles.text13,'String','Saved!','BackgroundColor','green')
                        set(handles.uipanel1,'BackgroundColor','green')
                     end
                  end
                  fclose(fid_a);
                  
               catch
                  errordlg('There is no log file. Please create "linescanlog.txt".',...
                     'linescanlog.txt Not Found!','modal')
               end
               
               fid_s = fopen([name,'.txt'],'r'); %% check to see if data previously saved
               if fid_s == -1
                  set(handles.text17,'String','Not Saved','BackgroundColor','red')
                  set(handles.uipanel3,'BackgroundColor','red')
               elseif fid_s >= 3
                  set(handles.text17,'String','Saved!','BackgroundColor','green')
                  set(handles.uipanel3,'BackgroundColor','green')
                  fclose(fid_s);
               end
               
               if n_out == 0
                  set(handles.text13,'String','Not Saved','BackgroundColor','red')
                  set(handles.uipanel1,'BackgroundColor','red')
                  
                  x = [handles.X(1),handles.X(handles.vectsize)]; %% approximate center of part
                  handles.XC = round(10*mean(x))/10;
                  y = [handles.Y(1),handles.Y(handles.vectsize)];
                  handles.YC = round(10*mean(y))/10;
               end
               
               maxx = handles.XC + 5.00;
               minx = handles.XC - 5.00;
               maxy = handles.YC + 5.00;
               miny = handles.YC - 5.00;
               
               set(handles.text7,'String',num2str(minx));
               set(handles.slider2,'Min',minx,'Max',maxx,'Value',handles.XC)
               set(handles.text8,'String',num2str(maxx));
               set(handles.text6,'String',num2str(handles.XC));
               
               set(handles.text10,'String',num2str(miny));
               set(handles.slider1,'Min',miny,'Max',maxy,'Value',handles.YC)
               set(handles.text11,'String',num2str(maxy));
               set(handles.text9,'String',num2str(handles.YC));
               
               handles.plotfigure = figure('Position',[6 65 1271 270],'Color','w');
               
               handles.dR = str2num(get(handles.edit4,'String'));
               plot_data(handles)
         end
      end
   end
end
guidata(handles.figure1,handles);
 
%%% Populates Listbox With Desired Filetypes
function load_listbox(dir_path,handles)
cd (dir_path); %% changes working directory path to direction in which program was opened
dir_struct = dir([dir_path '\*.']); %% displays "up one directory" symbols
dir_struct2 = dir([dir_path '\*.DAT']); %% displays all files with .DAT extention
dir_struct = [dir_struct; dir_struct2]; 
[sorted_names,sorted_index] = sortrows({dir_struct.name}'); %% sort files
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = [sorted_index];
guidata(handles.figure1,handles);
set(handles.listbox1,'String',handles.file_names,'Value',1)
set(handles.text1,'String',pwd)
 
%%% Plots Linescan Data Using Nominal Curve as Reference
function plot_data(handles)
R = -1364.75;  %% vertex radius in mm -1364.75 for modr
K = -2.8889;   %% conic constant
D = 840;      %% diameter in mm
A2 = +1.21993e-16;  %% aspheric rho^6 term
 
if get(handles.radiobutton2,'Value')
   R = R - handles.probedia/2;  %% correction for probe for pick-point data
end
 
R = R + handles.dR;
c = 1/R;
 
try
   for l = 1:handles.vectsize %% translates X and Y into radial positions
      if handles.Y(l) < handles.YC
         r(l,1) = -((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      else
         r(l,1) = ((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      end
   end
 
   Zc = (A2*r.^6) + (c*r.^2)./(1 + (1 - (K + 1)*c^2*r.^2).^0.5); %% calculate nominal surface
   
   diff = detrend(handles.Z - Zc); %% remove tilt
   
   midpoint = round(median(find((-6 < r & r < 6),1))); %% shift data down in the center
   piston = diff(midpoint);
   diff = diff - piston + eps;
   
catch
   handles.XC = 322.8; 
   handles.YC = 388.1;
   
   for l = 1:handles.vectsize %% translates X and Y into radial positions
      if handles.Y(l) < handles.YC
         r(l,1) = -((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      else
         r(l,1) = ((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      end
   end
   
   Zc = (A2*r.^6) + (c*r.^2)./(1 + (1 - (K + 1)*c^2*r.^2).^0.5); %% calculate nominal surface
   
   diff = detrend(handles.Z - Zc); %% remove tilt
   
   midpoint = round(median(find((-1.5 < r & r < 1.5),1))); %% shift data down in the center
   piston = diff(midpoint);
   diff = diff - piston + eps;
end
 
switch handles.hflip %% flips data horizontally
   case 1
      r = -r;
end
 
diff = diff + handles.piston;
 
if get(handles.togglebutton2,'Value') == 1 %% creats max/min values
   liney3 = 1000*mean(diff(find(diff == max(diff(1:round(size(diff,1)/2),1)))-4:...
      find(diff == max(diff(1:round(size(diff,1)/2),1)))+4));
   liney4 = 1000*mean(diff(find(diff == max(diff(round(size(diff,1)/2):...
      size(diff,1),1)))-4:find(diff == max(diff(round(size(diff,1)/2):size(diff,1),1)))+4));
   liney5 = 1000*mean(diff(find(diff == min(diff(1:round(size(diff,1)/2),1)))-4:...
      find(diff == min(diff(1:round(size(diff,1)/2),1)))+4));
   liney6 = 1000*mean(diff(find(diff == min(diff(round(size(diff,1)/2):...
      size(diff,1),1)))-4:find(diff == min(diff(round(size(diff,1)/2):size(diff,1),1)))+4));
end
 
if get(handles.radiobutton3,'Value') == 1 %% sets plot units
   linex = -500:1:500;
elseif get(handles.radiobutton4,'Value') == 1
   linex = -500/25.4:1:500/25.4;
   r = r/25.4;
elseif get(handles.radiobutton5,'Value') == 1
   linex = -500/415:.1:500/415;
   r = r/415;
end
 
if str2num(get(handles.edit9,'String')) == 0 %% smoothing and slope function
else
   smooth = str2num(get(handles.edit9,'String')); %% odd only
   for pp = 1:(size(r,1)-smooth+1)
      newdiff(pp,1) = (sum(diff(pp:(pp+smooth-1))))/smooth;
   end
   diff = [];
   for ppp = 1:smooth
      if ppp < (smooth/2 + .5)
         diff = [diff;newdiff(1)];
      elseif ppp == (smooth/2 + .5)
         diff = [diff;newdiff];
      elseif ppp > (smooth/2 + .5)
         diff = [diff;newdiff(size(newdiff,1))];
      end
   end
end
 
if str2num(get(handles.edit10,'String')) == 0 %% plot data without high slopes
else
   qq = 1;
   slopesize = str2num(get(handles.edit10,'String'));
   slopediff(1,1) = NaN;
   for i=1:slopesize:(size(r,1)-slopesize)
      qq = qq + 1;
      slopediff(qq,1) = (diff(i+slopesize)-diff(i))/(r(i+slopesize)-r(i));
      if abs((diff(i+slopesize)-diff(i))/(r(i+slopesize)-r(i))) >= 4.0E-4
         for ii = 0:slopesize-1
            diff(i+ii) = NaN;
         end
      end
   end
end
 
figure(1) %% plot
close(figure(2))
set(gcf,'Color','w')
axes('DrawMode','fast','YGrid','on','FontWeight','bold')
axis([min(linex) max(linex) str2num(get(handles.edit8,'String')) ...
   str2num(get(handles.edit7,'String'))])
xlabel('Radial Position [mm]')
ylabel('Surface Error [um]')
box('on')
hold('all')
 
switch get(handles.togglebutton6,'Value') %% plot slopes
   case 1
      plot(r,1000*diff,'b-',r(1:slopesize:size(r,1)),10000*slopediff,'r')
   case 0
      switch get(handles.togglebutton2,'Value') %% plot min/max lines
         case 1
            switch get(handles.togglebutton3,'Value') %% plot opd surface data
               case 1
                  switch get(handles.togglebutton4,'Value') %% plot earlier line
                     case 1
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxx,...
                                 handles.yyy,'m-',handles.xxxx,handles.yyyy,'k-',linex3,...
                                 liney,'r',linex4,liney,'r',linex5,liney,'r',linex6,liney,'r');
                           case 0
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxx,...
                                 handles.yyy,'m-',linex3,liney,'r',linex4,liney,'r',linex5,...
                                 liney,'r',linex6,liney,'r');
                        end
                     case 0
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxxx,...
                                 handles.yyyy,'k-',linex3,liney,'r',linex4,liney,'r',linex5,...
                                 liney,'r',linex6,liney,'r');
                           case 0
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',linex3,...
                                 liney,'r',linex4,liney,'r',linex5,liney,'r',linex6,liney,'r');
                        end
                  end
                  
               case 0
                  switch get(handles.togglebutton4,'Value') %% plot earlier line
                     case 1
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xxx,handles.yyy,'m-',...
                                 handles.xxxx,handles.yyyy,'k-',linex3,liney,'r',linex4,...
                                 liney,'r',linex5,liney,'r',linex6,liney,'r');
                           case 0
                              plot(r,1000*diff,'b-',handles.xxx,handles.yyy,'m-',...
                                 linex3,liney,'r',linex4,liney,'r',linex5,liney,'r',linex6,...
                                 liney,'r');
                        end
                     case 0
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xxxx,handles.yyyy,'k-',linex3,...
                                 liney,'r',linex4,liney,'r',linex5,liney,'r',linex6,liney,'r');
                           case 0
                              plot(r,1000*diff,'b-',linex,liney3,'r',linex,liney4,'r',linex,...
                                 liney5,'r',linex,liney6,'r');
                        end
                  end
            end
            
         case 0
            switch get(handles.togglebutton3,'Value') %% plot opd surface data
               case 1
                  switch get(handles.togglebutton4,'Value') %% plot earlier line
                     case 1
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxx,...
                                 handles.yyy,'m-',handles.xxxx,handles.yyyy,'k-');
                           case 0
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxx,...
                                 handles.yyy,'m-');
                        end
                     case 0
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r-',handles.xxxx,...
                                 handles.yyyy,'k-');
                           case 0
                              plot(r,1000*diff,'b-',handles.xx,handles.yy,'r.-');
                        end
                  end
               case 0
                  switch get(handles.togglebutton4,'Value') %% plot earlier line
                     case 1
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xxx,handles.yyy,'m-',...
                                 handles.xxxx,handles.yyyy,'k-');
                           case 0
                              plot(r,1000*diff,'b-',handles.xxx,handles.yyy,'m-');
                        end
                     case 0
                        switch get(handles.togglebutton5,'Value') %% plot piston alignment
                           case 1
                              plot(r,1000*diff,'b-',handles.xxxx,handles.yyyy,'k-');
                           case 0
                              plot(r,1000*diff,'b-');
                        end
                  end
            end
      end
end
guidata(handles.figure1,handles);
 
%%% Creates Listbox for File Selection
function listbox1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
   set(hObject,'BackgroundColor','white');
else
   set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
 
%%% Changes Y-Center of Linescan Data to Align Out Coma Using Slider
function slider1_Callback(hObject, eventdata, handles)
set(handles.text9,'String',num2str(get(handles.slider1,'Value')));
handles.YC = get(handles.slider1,'Value'); %% in units of mm
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates Y-Center Slider
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
 
%%% Changes X-Center of Linescan Data to Align Out Coma Using Slider
function slider2_Callback(hObject, eventdata, handles)
set(handles.text6,'String',num2str(get(handles.slider2,'Value')));
handles.XC = get(handles.slider2,'Value'); %% in units of mm
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates X-Center Slider
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
 
%%% Gets Value Entered into Angle Textbox
function edit1_Callback(hObject, eventdata, handles)
handles.angle = str2num(get(handles.edit1,'String'));
guidata(handles.figure1,handles);
 
%%% Creates Angle Textbox
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Horizontally Flips Linescan Data
function togglebutton1_Callback(hObject, eventdata, handles)
handles.hflip = get(hObject,'Value');
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Save Alignment Settings to a Log for Future Program Access
function pushbutton1_Callback(hObject, eventdata, handles)
if handles.angle < 0
   errordlg('There is no angle. Please input an angle.','Enter an Angle First','modal')
else
   fid_a = fopen('C:\MATLAB_SV701\work\MODS\linescanlog.txt','a+');
   fprintf(fid_a, '%11s,%14s,%6.3f,%6.3f,%1.0f,%3.0f,%1.0f\n',date,handles.file1,...
      handles.XC,handles.YC,handles.hflip,handles.angle,1); %% saves all alignment data
   fclose(fid_a);
   
   
   set(handles.text13,'String','Saved!','BackgroundColor','green')
   set(handles.uipanel1,'BackgroundColor','green') %% set text and color to saved setting
end
guidata(handles.figure1,handles);
 
%%% Saves Aligned Linescan Data to a File
function pushbutton2_Callback(hObject, eventdata, handles)
if handles.angle < 0
   errordlg('There is no angle. Please input an angle.','Enter an Angle First','modal')
else
   R = -1364.75;  %% vertex radius in mm -1364.75 for modr
   K = -2.8889;   %% conic constant
   D = 840;      %% diameter in mm
   A2 = +1.21993e-16;  %% aspheric rho^6 term
   if get(handles.radiobutton2,'Value')
      R = R - handles.probedia/2;  %% correction for probe for pick-point data
   end
   
   R = R + handles.dR;
   c = 1/R;
 
   for l = 1:handles.vectsize %% translates X and Y into radial positions
      if handles.Y(l) < handles.YC
         r(l,1) = -((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      else
         r(l,1) = ((handles.X(l) - handles.XC).^2 + (handles.Y(l) - handles.YC).^2).^0.5;
      end
   end
   
   Zc = (A2*r.^6) + (c*r.^2)./(1 + (1 - (K + 1)*c^2*r.^2).^0.5); %% calculate nominal surface
   
   diff = detrend(handles.Z - Zc); %% remove tilt
   
   midpoint = round(median(find((-6 < r & r < 6),1))); %% shift data down in the center
   piston = diff(midpoint);
   diff = diff - piston + eps;
   diff = diff + handles.piston;
   
   switch handles.hflip %% flips data horiztonally
      case 1
         r = -r;
   end
   
   if str2num(get(handles.edit9,'String')) == 0 %% smoothing and slope function
   else
      smooth = str2num(get(handles.edit9,'String')); %% odd only
      for pp = 1:(size(r,1)-smooth+1)
         newdiff(pp,1) = (sum(diff(pp:(pp+smooth-1))))/smooth;
      end
      diff = [];
      for ppp = 1:smooth
         if ppp < (smooth/2 + .5)
            diff = [diff;newdiff(1)];
         elseif ppp == (smooth/2 + .5)
            diff = [diff;newdiff];
         elseif ppp > (smooth/2 + .5)
            diff = [diff;newdiff(size(newdiff,1))];
         end
      end
   end
   
   if str2num(get(handles.edit10,'String')) == 0 %% plot data without high slopes
   else
      qq = 1;
      slopesize = str2num(get(handles.edit9,'String'));
      slopediff(1,1) = NaN;
      for i=1:slopesize:(size(r,1)-slopesize)
         qq = qq + 1;
         handles.slopediff(qq,1) = (diff(i+slopesize)-diff(i))/(r(i+slopesize)-r(i));
         if abs((diff(i+slopesize)-diff(i))/(r(i+slopesize)-r(i))) >= 4.0E-4 %% units in mrad
            for ii = 0:slopesize-1
               diff(i+ii) = NaN;
            end
         end
      end
   end
   
   theta=pi()*handles.angle/180; %% convert into output data
   data(:,1)=(r*sin(theta))/1000;
   data(:,2)=(r*cos(theta))/1000;
   data(:,3)=(diff*1000)/1000000;
   
   rdata(:,1)=r - min(r);
   rdata(:,2)=(diff*1000);
   
   switch get(handles.checkbox3,'Value') %% remove inner 30 mm from data
      case 0
         cutdata = 0;
      case 1
         cutdata = 30;
   end
   
   m=1;
   n=1;
   for i=1:size(r,1)
      if m<size(r,1)
         if r(m)<-cutdata
            data2(n,:) = data(m,:);
            n = n + 1;
         end
         if r(m)>cutdata
            data2(n,:) = data(m,:);
            n = n + 1;
         end
      end
      
      if get(handles.radiobutton1,'Value')
         m = m + 1;
      else
         m = m + 1;
      end
   end
   
   file3=[handles.file1,'.txt']; %% write out file
   file4=[handles.file1 '_r' '.txt'];
   dlmwrite(file3,data2,'-append','delimiter','\t');
   dlmwrite(file4,rdata,'-append','delimiter','\t');
   
   bigdatafile = ['MODSred_',handles.file1(1:6),'_AllData','.txt'];
   dlmwrite(bigdatafile,data2,'-append','delimiter',',');
   
   set(handles.text17,'String','Saved!','BackgroundColor','green')
   set(handles.uipanel3,'BackgroundColor','green') %% set text and color to saved setting
end
guidata(handles.figure1,handles);
 
%%% Display or Remove Max/Min Value Lines in Plot
function togglebutton2_Callback(hObject, eventdata, handles)
if get(handles.togglebutton2,'Value') == 1
   set(handles.togglebutton2,'String','Remove Max/Min Lines in Plot');
   plot_data(handles)
else
   set(handles.togglebutton2,'String','Plot Max/Min Lines');
   plot_data(handles)
end
guidata(handles.figure1,handles);
 
%%% Display or Remove Optical Surface Map Line Profile Data in Plot
function togglebutton3_Callback(hObject, eventdata, handles)
if get(handles.togglebutton3,'Value') == 1
   set(handles.togglebutton3,'String','Remove OPD Surface Data');
   
   fileOPD = [handles.file1,'_OPD.txt']; %% import optical profile data from PhaseMOSAIC
   [xx, yy] = textread(fileOPD,'%f %f','delimiter','\t');
   datasize = max(xx) - min(xx); %% gets size of data
   xx = xx - datasize/2; %% centers data about zero
   
   midpoint = round(median(find((-2 < xx & xx < 2),1))); %% shift data down in the center
   piston = yy(midpoint);
   yy = yy - piston + 6.000001;
   
   handles.xx = xx; %% globalizes vector
   handles.yy = yy; %% globalizes vector
   
   plot_data(handles)
   guidata(handles.figure1,handles);
else
   set(handles.togglebutton3,'String','Plot OPD Surface Data');
   plot_data(handles)
end
guidata(handles.figure1,handles);
 
%%% Closes Main Figure Window
function figure1_DeleteFcn(hObject, eventdata, handles)
close(figure(1))
 
%%% Indicates Whether to Remove Inner 157mm From Saved Data
function checkbox1_Callback(hObject, eventdata, handles)
 
%%% Pistons Optical Profile Data Up
function pushbutton3_Callback(hObject, eventdata, handles)
handles.yy = handles.yy + 0.2; %% in units of um
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Pistons Earlier Line Up
function pushbutton4_Callback(hObject, eventdata, handles)
handles.yyy = handles.yyy + 0.5; %% in units of um
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Pistons Earlier Line Down
function pushbutton5_Callback(hObject, eventdata, handles)
handles.yyy = handles.yyy - 0.5; %% in units of um
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Pistons Optical Profile Data Down
function pushbutton6_Callback(hObject, eventdata, handles)
handles.yy = handles.yy - 0.2; %% in units of um
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Text Input to Change Power of Linescan
function edit4_Callback(hObject, eventdata, handles)
handles.dR = str2num(get(handles.edit4,'String'));
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates Power Textbox
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Text Input to Indicate Which Earlier Date to Plot a Linescan of at the Same Angle
function edit5_Callback(hObject, eventdata, handles)
 
%%% Creates Earlier Date Textbox
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Plots Linescan of an Earlier Date with the Same Angle
function togglebutton4_Callback(hObject, eventdata, handles)
if get(handles.togglebutton4,'Value') == 1
   set(handles.togglebutton4,'String','Remove Earlier Line');
   edate = get(handles.edit5,'String');
   numeros = size(handles.file1,2); %% retrieves size of numerical filename
   orient = handles.file1(numeros-2:numeros); %% retrieves angular orientation of file
   
   try %% attempts to find one version of the previous file, if no find then finds another
      fileR = [edate,'_',orient,'_rv3.txt'];
      tet = 1;
      [xx, yy] = textread(fileR,'%f %f','delimiter','\t');
   catch
      fileR = [edate,'_',orient,'_r.txt'];
      tet = 0;
      [xx, yy] = textread(fileR,'%f %f','delimiter','\t');
   end
   
   datasize = max(xx) - min(xx); %% determines size of data set
   xx = xx - datasize/2; %% centers data about zero
   
   if tet == 0 %% shift data down in the center if particular data type was found
      midpoint = round(median(find((-2 < xx & xx < 2),1)));
      piston = yy(midpoint);
      yy = yy - piston + 0.000001;
   end
      
   if get(handles.radiobutton3,'Value') == 1
      handles.xxx = xx;
   elseif get(handles.radiobutton4,'Value') == 1
      handles.xxx = xx/25.4;
   elseif get(handles.radiobutton5,'Value') == 1
      handles.xxx = xx/415;
   end
   
   handles.yyy = yy; %% globalizes vector
   
   plot_data(handles)
   guidata(handles.figure1,handles);
else
   set(handles.togglebutton4,'String','Plot Earlier Line');
   plot_data(handles)
end
guidata(handles.figure1,handles);
 
%%% Indicates Whether to Save Linescan Data Without Inner 30mm
function checkbox3_Callback(hObject, eventdata, handles)
 
%%% Text Input to Choose a Filename to Plot Over Current Linescan
function edit6_Callback(hObject, eventdata, handles)
 
%%% Creates Filename Textbox
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Plots File Indicated by Textbox (edit6) Over Current Linescan for Alignment Purposes
function togglebutton5_Callback(hObject, eventdata, handles)
if get(handles.togglebutton5,'Value') == 1
   set(handles.togglebutton5,'String','Remove Earlier Line');
   fileR = get(handles.edit6,'String');
   [xx, yy] = textread([fileR '_r.txt'],'%f %f','delimiter','\t'); %% get previous data
   datasize = max(xx) - min(xx); %% gets data lateral size
   xx = xx - datasize/2; %% center data about zero
   for i = 1:size(yy,1)
      yy(i) = yy(i)-(.00/10^5)*(7*(xx(i))^2 - 1); %% change power of previous data if needed
   end
 
   midpoint = round(median(find((-2 < xx & xx < 2),1)));
   piston = yy(midpoint);
   yy = yy - piston + eps; %% shift data down in the center
   
   handles.xxxx = xx; %% globalizes vector
   handles.yyyy = yy; %% globalizes vector
   
   plot_data(handles)
   guidata(handles.figure1,handles);
else
   set(handles.togglebutton5,'String','Plot Earlier Line');
   plot_data(handles)
end
guidata(handles.figure1,handles);
 
%%% Pistons Current Linescan Data Up
function pushbutton8_Callback(hObject, eventdata, handles)
handles.piston = handles.piston + 0.0001; %% in units of mm
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Pistons Current Linescan Data Down
function pushbutton9_Callback(hObject, eventdata, handles)
handles.piston = handles.piston - 0.0001; %% in units of mm
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Changes Maximum Vertical Plot Scale
function edit7_Callback(hObject, eventdata, handles)
close(figure(1))
handles.plotfigure = figure('Position',[6 65 1271 270]);
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates Plot Scale Max Textbox
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Changes Minimum Vertical Plot Scale
function edit8_Callback(hObject, eventdata, handles)
close(figure(1))
handles.plotfigure = figure('Position',[6 65 1271 270]);
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates Plot Scale Min Textbox
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Indicates Window Size of a Moving Average to Smooth Linescan Data (odd number only)
function edit9_Callback(hObject, eventdata, handles)
if mod(str2num(get(handles.edit9,'String')),2) == 1 | str2num(get(handles.edit9,'String')) == 0
   plot_data(handles)
else
   errordlg('This number must be odd.','Enter an Odd Number','modal')
end
guidata(handles.figure1,handles);
 
%%% Creates Smoothing Window Size Textbox
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 
%%% Indicates "data point separation" Size for Which To Evaluate the Slope Over
function edit10_Callback(hObject, eventdata, handles)
plot_data(handles)
guidata(handles.figure1,handles);
 
%%% Creates Slope Textbox
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
 	
%%% Plot Calculated Slopes of Data Over Current Linescan
function togglebutton6_Callback(hObject, eventdata, handles)
if get(handles.togglebutton6,'Value') == 1
   set(handles.togglebutton6,'String','Remove Slopes');
   plot_data(handles)
   guidata(handles.figure1,handles);
else
   set(handles.togglebutton6,'String','Plot Slopes');
   plot_data(handles)
   guidata(handles.figure1,handles);
end
