function varargout = peregrine(varargin)
% PEREGRINE MATLAB code for peregrine.fig
%      PEREGRINE, by itself, creates a new PEREGRINE or raises the existing
%      singleton*.
%
%      H = PEREGRINE returns the handle to a new PEREGRINE or the handle to
%      the existing singleton*.
%
%      PEREGRINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEREGRINE.M with the given input arguments.
%
%      PEREGRINE('Property','Value',...) creates a new PEREGRINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before peregrine_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to peregrine_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help peregrine

% Last Modified by GUIDE v2.5 13-Jan-2014 17:00:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @peregrine_OpeningFcn, ...
                   'gui_OutputFcn',  @peregrine_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before peregrine is made visible.
function peregrine_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to peregrine (see VARARGIN)

% Choose default command line output for peregrine
handles.output = hObject;
handles.offline=0;

global debug
global timeit
global k
global kp
global ss 
global live1
live1=1;
set(handles.live, 'value', 1);
global frt
frt=0;
kp=1;


%%%%%%%%%%% Debug or performance tuning %%%%%%%%%%%%
debug=0;
timeit=0;
%%%%%%%%%%% Debug or performance tuning %%%%%%%%%%%%

ss=0;
k=1;
global view_type
view_type=1;
set(handles.record, 'enable', 'off');
if size(varargin,2)
    handles.offline=1; 
    handles.frmloc=[varargin{1}, '/'];
    handles=get_ready(handles, []);
end


global rp_id
global zoom_val
rp_id=0;
zoom_val=0;

instr1=sprintf('[I]\n(1)Use file > open to select a frame (or video)\n');
instr1=sprintf('%s(2)Press start to view the frames', instr1);
set_instr(instr1, handles);

%%%%%% keyboard shortcuts
set(handles.figure1, 'KeyPressFcn', @keypress_nav);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes peregrine wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = peregrine_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function bin1_Callback(hObject, eventdata, handles)
% hObject    handle to bin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bin_t=round(get(hObject,'Value'));
set(handles.bin1_txt, 'string', sprintf('%.3d',bin_t));


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function bin1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function area1_Callback(hObject, eventdata, handles)
% hObject    handle to area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
area1_t=round(get(hObject,'Value'));
set(handles.area1_txt, 'string', sprintf('%.3d', area1_t));

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function area1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to area1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function img1_Callback(hObject, eventdata, handles)
% hObject    handle to img1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img1_t=round(get(hObject,'Value'));
set(handles.img1_txt, 'string', sprintf('%.3d',img1_t));

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function img1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in bgyes.
function bgyes_Callback(hObject, eventdata, handles)
% hObject    handle to bgyes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgyes


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in handles.startstop.
function startstop_Callback(hObject, eventdata, handles)
% hObject    handle to handles.startstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ss
ss=~ss;

if ss 
    if handles.offline==3
        start(handles.vid);
    end
    set(handles.startstop, 'string', 'STOP');
elseif ~ss 
   if handles.offline==3
       flushdata(handles.vid);
       stop(handles.vid);
   end
   set(handles.startstop, 'string', 'START');
end

handles=pgn_loopdeloop(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in roi.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global roi_crop;
roi_crop=ceil(getrect(handles.axes1));

guidata(hObject, handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global roi_crop
global roi_cut

handles.conf(2,1)=round(get(handles.img1, 'value'));
handles.conf(2,2)=round(get(handles.area1, 'value'));
handles.conf(2,3)=round(get(handles.bin1, 'value'));
handles.conf(3,:)=roi_crop;
handles.conf(4:13,:)=roi_cut;

csvwrite([handles.frmloc, 'config_', handles.suff,  '.csv'], handles.conf);
if sum(handles.datfile(:,1))
    csvwrite([handles.frmloc, 'X00_', handles.suff, '.csv'], handles.datfile);
    datstr.frame=handles.frame;
    datstr.conf=handles.conf;
    save([handles.frmloc, 'dat0_', handles.suff, '.mat'], 'datstr');
end
set_instr('[I] data files saved...', handles);

% --- Executes on button press in cut.
function cut_Callback(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global roi_cut
nrz=find(roi_cut(:,1)==0);
if ~isempty(nrz)
    roi_cut(nrz(1),:)=ceil(getrect(handles.axes1));
else
    set(handles.instr, 'string', '[!] no more space available.');
end


% --- Executes when selected object is changed in viewtype.
function viewtype_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in viewtype 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global view_type
global kp    

set(handles.speed_up, 'enable', 'on');
set(handles.slow_down, 'enable', 'on');

if get(handles.view, 'value')
    view_type=1;
    set(handles.instr, 'string', '');
    set(handles.record, 'value', 0);
    set(handles.record, 'enable', 'off');
elseif get(handles.Mark, 'value')
    view_type=2;
    set(handles.instr, 'string', '');
    set(handles.record, 'value', 0);
    set(handles.record, 'enable', 'off');
    set(handles.status, 'foregroundcolor', 'r');
    if handles.shape
        set(handles.status, 'tooltip', ...
            'Number of targets and shape computation performance. :( means bad, fix thresholds, check background');
    else
        set(handles.status, 'tooltip', ...
            'Number of targets. Keep it consistent.');
    end
    set_instr('[I] Ensure that the number of targets are consistent throughout the video', handles);
elseif get(handles.track, 'value')
    view_type=3;
    set(handles.instr, 'string', '');
    set(handles.record, 'enable', 'on');
    if isfield('handles', 'sides_cm')
        if ~handles.sides_cm(1)
            set_instr('[!] Region not calibrated in preferences', handles);
        else
            set_instr('[I] Check the box if you wish to record', handles);
        end
    end
    % don't change speeds during tracking ....
    kp=1;
    set(handles.speed_up, 'enable', 'off');
    set(handles.slow_down, 'enable', 'off');
    set(handles.status, 'foregroundcolor', 'g');
    set(handles.status, 'tooltip', ...
        'Percentage of continuous tracks (should be high) /Number of tracks (should be low)');
elseif get(handles.repair, 'value')
    view_type=4;
    set(handles.record, 'value', 0);
    set(handles.record, 'enable', 'off');
    set(handles.instr, 'string', '');
    set(handles.status, 'foregroundcolor', 'b');
    set(handles.status, 'tooltip', 'target id (% done)');  
elseif get(handles.verify, 'value')
    view_type=5;
    set_instr('Only recorded tracks are shown. Images are saved in /tmp', handles);
    set(handles.record, 'value', 0);
    set(handles.record, 'enable', 'off');
end


% --------------------------------------------------------------------
function m_file_Callback(hObject, eventdata, handles)
% hObject    handle to m_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_about_Callback(hObject, eventdata, handles)
% hObject    handle to m_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aboutstr=[];

aboutstr=[aboutstr sprintf('Peregrine (v 0.1) ... \n')];
aboutstr=[aboutstr sprintf('... is a falcon that hunts flocking birds\n')];
aboutstr=[aboutstr sprintf('... is a multi-target tracker built using MATLAB\n\n')];

aboutstr=[aboutstr sprintf('Sachit Butail, January 13, 2013, New York, USA\n')];
aboutstr=[aboutstr sprintf('================================\n')];

aboutstr=[aboutstr sprintf('***Credits***:\n')];
aboutstr=[aboutstr sprintf('- Dynamical Systems Laboratory, Polytechnic Institute of NYU\n')];
aboutstr=[aboutstr sprintf('- Collective Dynamics and Control Laboratory, University of Maryland\n')];
aboutstr=[aboutstr sprintf('- The following files are downloaded from MATLAB Central File Exchange: \n')];
aboutstr=[aboutstr sprintf('munkres.m \ninputsdlg.m \nemgm.m \n')];

msgbox(aboutstr, 'About');

% --------------------------------------------------------------------
function m_edit_Callback(hObject, eventdata, handles)
% hObject    handle to m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=change_prefs(handles);
guidata(hObject, handles);
% --------------------------------------------------------------------
function m_help_Callback(hObject, eventdata, handles)
% hObject    handle to m_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --------------------------------------------------------------------
function m_open_Callback(hObject, eventdata, handles)
% hObject    handle to m_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k
k=1;
global ss 
ss=0;
global view_type
view_type=1;

[filename, pathname, filterindex]=uigetfile(...
            {'*.bmp; *.tif; *.png; *.jpg; *.jpeg;', 'Image files'; ...
             '*.avi; *.mov;', 'Video files'}, 'Choose an image file or video');
if filename

    handles.offline=filterindex;
    handles.frmloc=[pathname, '/'];
    
    handles=get_ready(handles, filename);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function m_camera_Callback(hObject, eventdata, handles)
% hObject    handle to m_camera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
global k
k=1;
global ss 
ss=0;
global view_type
view_type=1;

handles.offline=3;
handles.frmloc='./';

handles=get_ready(handles, []);

        
        
% Update handles structure
guidata(hObject, handles);


function handles = get_ready(handles, filename)
global roi_crop;
global roi_cut;
global k;
set(gcf, 'pointer', 'watch');
set(handles.instr, 'string', 'Loading data...');
[handles.conf, handles.vid, handles.getfrm, ...
        handles.nfrm handles.fps, handles.nt, ...
        handles.bg, handles.alpha, handles.smooth handles.blur, handles.bb_blur, handles.bb_size, ...
        handles.suff max_val imh imw handles.fgislight handles.trktype, ...
        handles.record_verify roi_crop roi_cut img_t, area_t, ...
        bin_t handles.sides_cm handles.circ, handles.shape, handles.split]= init_setup(handles.frmloc, handles.offline, filename);

    
set(handles.m_preferences, 'enable', 'on');   
set(handles.ui_edit_preferences, 'enable', 'on');
if isempty(handles.conf)
    return
end
%%% make the video clearer
%     src=getselectedsource(vid);
%     set(src, 'Saturation', 50);
%     set(src, 'Sharpness', 10);
%     set(src, 'Contrast', 3);
%     set(src, 'Brightness',150);

handles.X=[];
handles.P=[];
handles.frame(handles.nfrm,1).Zk=[];
handles.frame(handles.nfrm,1).Mk=[];
handles.frame(handles.nfrm,1).X=[];
handles.frame(handles.nfrm,1).P=[];

%%%%%%
% axes(handles.axes1);
% axis([1 imh 1 imw]);
% set(handles.axes1, 'units', 'pixels');
% set(handles.axes1, 'position', [5 80 imw imh]);

set(handles.img1, 'Min',1,'Max',max_val,...
      'SliderStep',[1/(max_val-1) ceil(max_val*.01)/(max_val-1)]);
set(handles.area1, 'Min',0,'Max',100,...
      'SliderStep',[1/(100) ceil(100*.01)/(100-1)]);      
set(handles.bin1, 'Min',0,'Max',max_val,...
          'SliderStep',[1/(max_val) ceil(max_val*.01)/(max_val-1)]); 
set(handles.frmpos, 'Min',1,'Max',handles.nfrm,...
          'SliderStep',[1/(handles.nfrm-1) ceil(handles.nfrm*.01)/(handles.nfrm-1)]);       

set(handles.img1, 'value', img_t);
set(handles.img1_txt, 'string', sprintf('%.3d', img_t));
set(handles.area1, 'value', area_t);
set(handles.area1_txt, 'string', sprintf('%.3d', area_t));
set(handles.bin1, 'value', bin_t);
set(handles.bin1_txt, 'string', sprintf('%.3d', bin_t));
set(handles.frmpos, 'Value',1);  

% if exist([handles.frmloc, 'Z01_', handles.suff, '.csv'], 'file')
%     handles.Z1=csvread([handles.frmloc, 'Z01_', handles.suff, '.csv']);
% else
%     handles.Z1=zeros(3*handles.nfrm, handles.nt);
% end

handles.sz1=zeros(1,1000);

mm1=floor(k/handles.fps/60);
ss1=floor((k-mm1*handles.fps*60)/handles.fps);
set(handles.time, 'string', sprintf('%.2d:%.2d (%.5d)', mm1, ss1, k));
    
%%%% file for saving tracks
csvfile=[handles.frmloc, 'X00_', handles.suff, '.csv'];
if exist(csvfile, 'file')
    handles.datfile=csvread(csvfile);
    
    %%%% Backward compatibility
    if size(handles.datfile,1) < handles.nfrm*handles.nt
        handles.datfile=[handles.datfile; ...
             ones(handles.nfrm*handles.nt-size(handles.datfile,1),1)*zeros(1,size(handles.datfile,2))];
    end
    set(handles.verify, 'enable', 'on');
    set(handles.repair, 'enable', 'on');
    set(handles.zoomin, 'enable', 'on');
    set(handles.zoomout, 'enable', 'on');
    instr1=sprintf('[I]Looks like you have tracked this dataset before. You may\n');
    instr1=sprintf('%s(1) Retrack with different thresholds\n', instr1);
    instr1=sprintf('%s(2) Repair or Verify', instr1);
    set_instr(instr1, handles);
else
    handles.datfile=zeros(handles.nfrm*handles.nt, 20);
    instr1=sprintf('[I]\n(1)Use a combination of thresholds and crop/cut to obtain the best background\n');
    instr1=sprintf('%s(2)Select Mark to see if the actual targets are marked\n', instr1);
    instr1=sprintf('%s(3)Pay attention to the number of targets on the status bar so that it is consistent\n', instr1);
    instr1=sprintf('%s(4)If satisfied, *SAVE* and move to the frame from where you wish to begin tracking, select Track and check the box next to it to record\n', instr1);
    set_instr(instr1, handles);
end

%%% new datafile
dat0file=[handles.frmloc, 'dat0_', handles.suff, '.mat'];
if exist(dat0file, 'file')
    tmpstr=load(dat0file);
    handles.frame=tmpstr.datstr.frame;
end
set(gcf, 'pointer', 'arrow');


% --- Executes on slider movement.
function frmpos_Callback(hObject, eventdata, handles)
% hObject    handle to frmpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global k
k=ceil(get(hObject, 'Value'));
handles=onestep_wrapper(handles); 
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

function handles=onestep_wrapper(handles)

global k
global view_type
global rp_id

img_t=get(handles.img1, 'value');
bin_t=get(handles.bin1, 'value');
area_t=get(handles.area1, 'value');
write_tracks=get(handles.record, 'value');
axes(handles.axes1);
[handles.bg handles.nZ handles.X handles.P handles.datfile handles.cm2pix handles.sz1 mesg handles.frame]=onestep(handles.sides_cm, img_t, bin_t, area_t, ...
                handles.fgislight, handles.circ, handles.smooth, handles.bg, handles.alpha,...
                handles.blur, handles.bb_blur, handles.bb_size, handles.X, handles.P, handles.fps, ...
                handles.datfile, handles.getfrm, get(handles.bgyes, 'value'), handles.trktype, ...
                write_tracks, handles.record_verify, handles.split, handles.shape, handles.sz1, handles.frame);

if view_type==2
    if isfield(handles, 'nZ')
        set(handles.status, 'string', sprintf('%d', handles.nZ));
    end
elseif view_type==4

    if ~rp_id
        set_instr('[I] Mark the fish you wish to verify/repair', handles);
    end
%     set(handles.status, 'string', sprintf('%d', rp_id));
    perc_repair_done=sum(handles.datfile(:,2)==rp_id)/...
            (max(handles.datfile(:,1))-min(handles.datfile(handles.datfile(:,1)~=0,1)));
    set(handles.status, 'string', sprintf('%d (%.1f%% done)', rp_id, perc_repair_done*100));
end            
mm1=floor(k/handles.fps/60);
ss1=floor((k-mm1*handles.fps*60)/handles.fps);
set(handles.time, 'string', sprintf('%.2d:%.2d (%.5d)', mm1, ss1, k));

% --- Executes during object creation, after setting all properties.
function frmpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frmpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in record.
function record_Callback(hObject, eventdata, handles)
% hObject    handle to record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k;
if get(hObject, 'Value')
    set_instr(sprintf('[I] Recording tracks into a file from frame %d', k), handles);
    set(handles.verify, 'enable', 'on');
end
% Hint: get(hObject,'Value') returns toggle state of record


% --- Executes on button press in live.
function live_Callback(hObject, eventdata, handles)
% hObject    handle to live (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of live
global live1
live1=get(hObject, 'Value');
if ~live1
    set(hObject, 'string', 'Offline');
else
    set(hObject, 'string', 'Live');
end


% --------------------------------------------------------------------
function copy_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to copy_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global roi_crop;
global roi_cut;


config_file=[handles.frmloc, 'config_', handles.suff, '.csv'];  

if exist(config_file, 'file')
    set_instr('[?] a config file already exists...', handles);
else
    [filename, pathname, filterindex]=uigetfile(...
            '*.csv;', 'Choose an existing preference file');
    if filename
        copyfile([pathname, '/', filename], config_file);
        handles.conf=csvread(config_file);
        
        handles.nt=handles.conf(1,1);
        handles.fps=handles.conf(1,3);

        handles.img_t=handles.conf(2,1);
        handles.area_t=handles.conf(2,2);
        handles.bin_t=handles.conf(2,3);

        roi_crop=handles.conf(3,:);
        roi_cut=handles.conf(4:13,:);  

        handles.smooth=handles.conf(14,:);
        handles.blur=handles.conf(15,:);
        handles.bb_blur=handles.conf(16,:);

        handles.fgislight=handles.conf(17,1);
        handles.trktype=handles.conf(17,2);
        handles.bb_size=handles.conf(17,3);
        handles.record_verify=handles.conf(17,4);

        handles.sides_cm=handles.conf(18,1:2); 
        handles.split=handles.conf(18,3);
        handles.circ=handles.conf(18,4);
        handles.shape=handles.conf(19,1);
        handles.alpha=handles.conf(19,2);
        
        
        set(handles.img1, 'value', handles.img_t);
        set(handles.img1_txt, 'string', sprintf('%.3d', handles.img_t));
        set(handles.area1, 'value', handles.area_t);
        set(handles.area1_txt, 'string', sprintf('%.3d', handles.area_t));
        set(handles.bin1, 'value', handles.bin_t);
        set(handles.bin1_txt, 'string', sprintf('%.3d', handles.bin_t));
        
        guidata(hObject, handles);
    end
end
    

function set_instr(str, handles)
% fprintf('[%d]%s\n', handles.k, str)

if strcmp(str(1:3), '[!]')
    set(handles.instr, 'BackgroundColor', [1 .25 .25]);
elseif strcmp(str(1:3), '[A]')
    set(handles.instr, 'BackgroundColor', 'b'); 
elseif strcmp(str(1:3), '[W]')
    set(handles.instr, 'BackgroundColor', [0.75 .25 .25 ]);
else
    set(handles.instr, 'BackgroundColor', [0 .75 0]);
end
if ~isfield(handles, 'suff')
    handles.suff=[];
end
str=sprintf('%s\n%s', handles.suff, str);
set(handles.instr, 'String', str, 'fontsize', 10);

% strold=get(handles.txt_instr1, 'String');
% 
% str=sprintf('%s [%d, %d] %s;', strold, handles.k, handles.mqid, str);
% 
% str_cell=textscan(str, '%s', 'Delimiter', ';');
% 
% set(handles.txt_instr1, 'String', sprintf('%s;', str_cell{1}{end-1:end}));


function handles=pgn_loopdeloop(handles)
global ss
global k
global view_type
global timeit
global kp

global rp_id

while ss
    if timeit, tic; end
    img_t=get(handles.img1, 'value');
    bin_t=get(handles.bin1, 'value');
    area_t=get(handles.area1, 'value');
    write_tracks=get(handles.record, 'value');
    if timeit, 
        fprintf('after getting values: %.2f\n', toc);
        tic
    end
    axes(handles.axes1);
    
    [handles.bg handles.nZ handles.X handles.P handles.datfile handles.cm2pix handles.sz1 mesg handles.frame]=onestep(handles.sides_cm, img_t, bin_t, area_t, ...
                    handles.fgislight, handles.circ, handles.smooth, handles.bg, handles.alpha,...
                    handles.blur, handles.bb_blur, handles.bb_size, handles.X, handles.P, handles.fps, ...
                    handles.datfile, handles.getfrm, get(handles.bgyes, 'value'), handles.trktype, ...
                    write_tracks, handles.record_verify, handles.split, handles.shape, handles.sz1, handles.frame);
    mm1=floor(k/handles.fps/60);
    ss1=floor((k-mm1*handles.fps*60)/handles.fps);
    set(handles.time, 'string', sprintf('%.2d:%.2d (%.5d)', mm1, ss1, k));
    
    
    if timeit, 
        fprintf('after onestep: %.2f\n', toc);
        tic
    end
    obs=handles.datfile(:,5);
    ksteps=handles.datfile(:,1);
    nt_k=numel(unique(handles.datfile(:,2)));
    if view_type==3
        set(handles.status, 'string', sprintf('%.1f%% / %d %s', sum(obs~=0)/sum(ksteps~=0)*100, nt_k, mesg.txt));
    elseif view_type==2
%         if isfield(handles, 'nZ')
%             set(handles.status, 'string', sprintf('%d %s', handles.nZ, mesg.txt));
%         end
        set(handles.status, 'string', sprintf('%d %s', size(handles.frame(k).Zk,2), mesg.txt));
        if size(handles.frame(k).Zk,2)~=handles.nt && get(handles.repair_posonly, 'value')
            ss=0;
            set(handles.startstop, 'string', 'START');
            instr1=sprintf('[!] Measurements are not equal to targets. Fix needed\n');
            instr1=sprintf('%s Possible fixes:\n', instr1);
            instr1=sprintf('%s(1) add a target (a)\n', instr1);
            instr1=sprintf('%s(2) delete a target (d)\n', instr1);
            set_instr(instr1, handles);
        end
    elseif view_type==4  
        if ~rp_id % if there is no repair id 
            set_instr('[I] Mark the fish you wish to verify/repair', handles);
            
        else
            perc_repair_done=sum(handles.datfile(:,2)==rp_id)/...
                (max(handles.datfile(:,1))-min(handles.datfile(handles.datfile(:,1)~=0,1)));
            set(handles.status, 'string', sprintf('%d (%.1f%% done)', rp_id, perc_repair_done*100));
        end
        pos=handles.datfile(handles.datfile(:,1)==k & handles.datfile(:,2)==rp_id,3:4);
        
        if isempty(pos) && rp_id
            kp=1;
            ss=0;
            set(handles.startstop, 'string', 'START');
            instr1=sprintf('[!] Target track ended at %d, Fix and restart.\n', k);
            instr1=sprintf('%s Possible fixes:\n', instr1);
            instr1=sprintf('%s(1) if the target is unmarked, add a point (a)\n', instr1);
            instr1=sprintf('%s(2) if the target is marked with a circle, update the track (s)\n', instr1);
%             instr1=sptrintf('%s 3) You may switch tracks after an occlusion that was not resolved correctly (s)', instr1);
            set_instr(instr1, handles);
        end
    end
   
    % keep looping unless you are tracking, in which case you should end
    % when the images end
    %     if handles.offline && ~get(handles.Mark, 'value')
    if handles.offline   
        if get(handles.record, 'Value') && k==handles.nfrm
            ss=0;
            set(handles.startstop, 'string', 'START');
        end
        if get(handles.verify, 'Value') && handles.record_verify && k==handles.nfrm
            ss=0;
            set(handles.startstop, 'string', 'START');
        end
        if k ==handles.nfrm, k=1; end
    end
    set(handles.frmpos, 'value', k);
    
    if ss % increment only if the ss is set to nonzero
        k=min(k+kp,handles.nfrm);
    end
    if timeit, 
        fprintf('after the rest: %.2f\n', toc);
    end
end

% --- Executes on button press in zoomin.
function zoomin_Callback(hObject, eventdata, handles)
% hObject    handle to zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zoom_val

zoom_val=min(zoom_val+1, 10);


% --- Executes on button press in zoomout.
function zoomout_Callback(hObject, eventdata, handles)
% hObject    handle to zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global zoom_val

zoom_val=max(zoom_val-1, 1);


%%%%%%% keyboard shortcuts
function keypress_nav(src, event)

% global k

% event.Key
% event.Modifier
% event.Character
handles = guidata(src);
if ~isempty(event.Modifier)
    if strcmp(event.Modifier{:}, 'control')
        switch event.Key
            case 's'
               save_Callback(handles.save, [], handles);
        end
    end
else
    switch event.Key
        case 't'
              select_rp_id_Callback(handles.select_rp_id, [], handles);
        case 'i'
              set_rp_id_Callback(handles.set_rp_id, [], handles);
        case 's'
              switch_tracks_Callback(handles.switch_tracks, [], handles);
        case 'd'
    %           nfdel=input('Frames in future to delete []=10, #=other: '); 
    %           if isempty(nfdel), nfdel=10; end
              nfdel=5;
              delete_track_Callback(handles.delete_track, nfdel, handles);
    %     case 'x'
    %           chop_t_Callback(handles.chop_t, [], handles);      
        case 'a'
              add_t_Callback(handles.add_t, [], handles);
        case 'm'
              measure_t_Callback(handles.measure_t, [], handles);              
        case 'space'
              startstop_Callback(handles.startstop, [], handles);
        case 'leftarrow'
              prev_frame_Callback(handles.prev_frame, [], handles);
        case 'rightarrow'
              next_frame_Callback(handles.next_frame, [], handles);
    end
end

% --- Executes on button press in switch_tracks.
function switch_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to switch_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rp_id
global k
global rp_handle
global ss
global roi_crop

val_t=ceil(max(roi_crop(3:4))/10);

set(handles.instr, 'string', '[I] Click on the target. (Esc. to Cancel)');

% click on second track (the one you want to switch with the current id)
[x y button]=ginput(1);
if button==1
    Xk=handles.datfile(handles.datfile(:,1)==k,:);
    pos=Xk(:,3:4);
    [xp yp]=handles.cm2pix(pos(:,1), pos(:,2));
    dist=sum(([x; y]*ones(1,numel(xp))-[xp'; yp']).^2);
    [val idx1]=min(dist);
    if ~isempty(val) && val<=val_t
        idx=Xk(idx1, 2);

        % switch ids
        sidx1=(handles.datfile(:,2)==rp_id & handles.datfile(:,1)>=k);
        sidx2=(handles.datfile(:,2)==idx & handles.datfile(:,1)>=k);

        handles.datfile(sidx1,2)=idx;
        handles.datfile(sidx2,2)=rp_id;

        if exist('rp_handle', 'var')
            if ishandle(rp_handle)
                delete(rp_handle);
            end
        end
        rp_handle=plot(xp(idx1), yp(idx1), 'k+');

        %%%%%%%%% update information for switching etc.
%         handles.datfile(handles.datfile(:,2)==rp_id, 20)=1;


        %%%% interpolate
        id1=find(handles.datfile(:,2)==rp_id & handles.datfile(:,1)<k & handles.datfile(:,3)~=0);
        [val idx]=max(handles.datfile(id1,1));
        id1=id1(idx);
        id2=find(handles.datfile(:,2)==rp_id & handles.datfile(:,1)==k);
        if numel(id2)>1, id2=id2(1);end

        k1=handles.datfile(id1,1); % last non-zero k
        k2=handles.datfile(id2,1); % current k

        %%%%%%%% interpolate between the two
        nfld=size(handles.datfile,2);
        curr=zeros(numel(k1:k2), nfld);
        for jj=1:nfld
            curr(:,jj)=interp1([k1 k2], handles.datfile([id1,id2], jj), k1:k2);
        end

        % delete all rows that are lie between
        del=handles.datfile(:,1)<=k2 & handles.datfile(:,1)>=k1 & handles.datfile(:,2)==rp_id;

        handles.datfile(del,:)=[];

        nzi=find(handles.datfile(:,1)~=0); % find non zero entries
        nzi=nzi(end);
        handles.datfile(nzi+1:nzi+size(curr,1),:)=curr;


        if ~ss
            set_instr(sprintf('[I] Monitor the target fish'), handles);
            ss=1;
            set(handles.startstop, 'string', 'STOP');
            handles=pgn_loopdeloop(handles);
        end

        guidata(hObject, handles);
    else
        instr1=sprintf('[!]The target that you selected is very far from the marked point (Try again)\n');
        instr1=sprintf('%s There is no marked target in this frame. Use (a) to add a point.', instr1);
        set_instr(instr1, handles);
    end
end

% --- Executes on button press in select_rp_id.
function select_rp_id_Callback(hObject, eventdata, handles)
% hObject    handle to select_rp_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rp_id
global k

[x y]=ginput(1);
Xk=handles.datfile(handles.datfile(:,1)==k,:);
pos=Xk(:,3:4);
[xp yp]=handles.cm2pix(pos(:,1), pos(:,2));
dist=sum(([x; y]*ones(1,numel(xp))-[xp'; yp']).^2);
[val idx]=min(dist);

%         if val < 10
rp_id=Xk(idx, 2);
plot(xp(idx), yp(idx), 'k+');
%         end


% --------------------------------------------------------------------
function ui_open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global k
k=1;
global ss 
ss=0;
global view_type
view_type=1;

[filename, pathname, filterindex]=uigetfile(...
            {'*.bmp; *.tif; *.png; *.jpg; *.jpeg;', 'Image files'; ...
             '*.avi; *.mov;', 'Video files'}, 'Choose an image file or video');
if filename

    handles.offline=filterindex;
    handles.frmloc=[pathname, '/'];
    
    handles=get_ready(handles, filename);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function ui_edit_preferences_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_edit_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=change_prefs(handles);
guidata(hObject, handles);

function handles=change_prefs(handles)

prompt={'Foreground is darker (0=no, 1=yes):';... % 1
        'Frame rate:'; ... % 2
        'Number of targets (approx.):'; ... %3
        'Blurring big blobs (width, height, sigma, extension):'; ... % 4
        'Blurring blob size (pixels):'; ... % 5
        'Tracker type (0=kf+gnn, 1=kf+gnn+nn, 2=pf+gnn) :'; ... % 6 
        'Record frames'; ... % 7 during verification
        'Region of Interest (ROI) size (length, breadth in cm)'; ... % 8
        'Shape of ROI (0=square, 1=circle):'; ... % 9 
        'Split occlusions (1=yes, 0=no):'; ... % 10 (slow and only if similar targets present)
        'Shape tracking (1=yes, 0=no. for fish only): '; ... %11 
        'Running background alpha (<1)'}; % build (any number less than 1=running, 0=highest intensity)
name='Preferences';

% -- begin
% formats(1,1) = struct('type','check','style', {'checkbox'}, 'items', {''}, 'format','integer','limits',[0 1]);
% formats(2,1) = struct('type','edit','style', 'edit','items', {''}, 'format','integer','limits',[1 100]);
% formats(3,1) = struct('type','edit','style', 'edit','items', {''}, 'format','integer', 'limits', [0 10]);
% formats(4,1) = struct('type','edit','style', 'edit','items', {''}, 'format','text','limits',[0 1]);
% formats(5,1) = struct('type','edit','style', 'edit','items', {''}, 'format','integer','limits',[0 1000]);
% formats(6,1) = struct('type','list','style', 'listbox','items', {{'kf+gnn'; 'kf+gnn+nn'; 'pf+gnn'}}, 'format','integer','limits',[0 3]);
% formats(7,1) = struct('type','check','style', {'checkbox'}, 'items', {''}, 'format','integer','limits',[0 1]);
% formats(8,1) = struct('type','edit','style', 'edit','items', {''}, 'format','text','limits',[0 1]);
% formats(9,1) = struct('type','list','style', 'radiobutton','items', {{'SQ'; 'CI'}}, 'format','integer','limits',[0 1]);
% formats(10,1) = struct('type','check','style', {'checkbox'}, 'items', {''}, 'format','integer','limits',[0 1]);
% formats(11,1) = struct('type','check','style', {'checkbox'}, 'items', {''}, 'format','integer','limits',[0 1]);
% formats(12,1) = struct('type','edit','style', 'edit','items', {''}, 'format','float','limits',[0 1]);
% for ii=1:size(formats,2)
%     formats(ii).size=0;
% end
% defaultanswer={ handles.fgislight; ...
%                 handles.fps; ...
%                 handles.nt; ...
%                 sprintf('%d,', handles.conf(16,:));...
%                 handles.bb_size; ...
%                 handles.trktype+1; ...
%                 handles.record_verify; ...
%                 sprintf('%d,', handles.sides_cm); ...
%                 handles.circ+1; ...
%                 handles.split; ...
%                 handles.shape; ...
%                 handles.alpha};
% % options.WindowStyle='modal';
% options.AlignControls='on';
% options.UnitsMargin=10;
% % options.Resize='on';
% answer=inputsdlg(prompt, name, formats, defaultanswer, options); 
% keyboard
% -- end

numlines=1;
if ~handles.bb_size, handles.bb_size=200; end
defaultanswer={ sprintf('%d', handles.fgislight), ...
                sprintf('%d', handles.fps), ...
                sprintf('%d', handles.nt), ...
                sprintf('%d,', handles.conf(16,:)),...
                sprintf('%d', handles.bb_size), ...
                sprintf('%d', handles.trktype), ...
                sprintf('%d', handles.record_verify), ...
                sprintf('%d,', handles.sides_cm), ...
                sprintf('%d', handles.circ), ...
                sprintf('%d', handles.split), ...
                sprintf('%d', handles.shape), ...
                sprintf('%.2f', handles.alpha)};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
    handles.fgislight=str2num(answer{1});
    handles.fps=str2num(answer{2});
    handles.nt=str2num(answer{3});
    handles.bb_blur=str2num(answer{4});
    handles.bb_size=str2num(answer{5});
    handles.trktype=str2num(answer{6});
    handles.record_verify=str2num(answer{7});
    handles.sides_cm=str2num(answer{8});
    handles.circ=str2num(answer{9});
    handles.split=str2num(answer{10});
    handles.shape=str2num(answer{11});
    handles.alpha=str2num(answer{12});
    
    handles.conf(17,1)=handles.fgislight;
    handles.conf(1,3)=handles.fps;
    handles.conf(1,1)=handles.nt;
    handles.conf(16,:)=handles.bb_blur;
    handles.conf(17,2)=handles.trktype;
    handles.conf(17,3)=handles.bb_size;
    handles.conf(17,4)=handles.record_verify;
    handles.conf(18,1:2)=handles.sides_cm;
    handles.conf(18,3)=handles.split;
    handles.conf(18,4)=handles.circ;
    handles.conf(19,1)=handles.shape;
    handles.conf(19,2)=handles.alpha;
end

% save the preferences
save_Callback(handles.save, [], handles);

if handles.shape && handles.trktype ~=2
    set_instr('[!] if tracking shape tracker type should be 2', handles);
end

% --- Executes on button press in no_zoom.
function no_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to no_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zoom_val
zoom_val=0;


% --- Executes on button press in set_rp_id.
function set_rp_id_Callback(hObject, eventdata, handles)
% hObject    handle to set_rp_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rp_id

prompt={'Set new fish id (current id shown)'};
name='Set new fish id';
numlines=1;

defaultanswer={sprintf('%d', rp_id)};
answer=inputdlg(prompt,name,numlines,defaultanswer);

new_id=str2num(answer{1});

% switch ids
sidx1=(handles.datfile(:,2)==rp_id );
sidx2=(handles.datfile(:,2)==new_id);

handles.datfile(sidx1,2)=new_id;
handles.datfile(sidx2,2)=rp_id;

rp_id=new_id;

guidata(hObject, handles);


% --- Executes on button press in slow_down.
function slow_down_Callback(hObject, eventdata, handles)
% hObject    handle to slow_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global kp

kp=max(1, kp-1);

% --- Executes on button press in speed_up.
function speed_up_Callback(hObject, eventdata, handles)
% hObject    handle to speed_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global kp
kp=min(handles.fps, kp+1);


% --- Executes on button press in delete_track.
function delete_track_Callback(hObject, eventdata, handles)
% hObject    handle to delete_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k
global roi_crop
global view_type
global ss

if view_type==2 && get(handles.repair_posonly, 'value')
    set(handles.instr, 'string', '[I] Click on the target. (Esc. to cancel)');

    % click on second track (the one you want to switch with the current id)
    [x y button]=ginput(1);
    if button==1   
        % if there is no marked positions then copy all from the current
        % set of measurements
        if isempty(handles.frame(k).Mk)
            handles.frame(k).Mk=handles.frame(k).Zk;
        end
        
        xp=handles.frame(k).Mk(1,:);
        yp=handles.frame(k).Mk(2,:);
        
        dist=sum(([x; y]*ones(1,numel(xp))-[xp; yp]).^2);
        [val idx]=min(dist);
        plot(xp(idx), yp(idx), 'ro', 'markersize', 10);
        
        handles.frame(k).Mk(:,idx)=[];
        
        % of expected targets same as measurements
        if ~ss && size(handles.frame(k).Mk,2)==handles.nt 
            ss=1;
            set(handles.startstop, 'string', 'STOP');
            handles=pgn_loopdeloop(handles);
        end
        
        guidata(hObject, handles);
    end
elseif view_type==4

    instr1=sprintf('[I] Click on the track to delete 5 frames from %d frame onwards. (Esc. to cancel)', k);
    set_instr(instr1, handles);

    val_t=ceil(max(roi_crop(3:4))/10);

    [x y button]=ginput(1);
    if button==1
        Xk=handles.datfile(handles.datfile(:,1)==k,:);
        pos=Xk(:,3:4);
        [xp yp]=handles.cm2pix(pos(:,1), pos(:,2));
        dist=sum(([x; y]*ones(1,numel(xp))-[xp'; yp']).^2);
        [val idx]=min(dist);
        if ~isempty(val) && val<=val_t

            del_id=Xk(idx, 2);
            plot(xp(idx), yp(idx), 'ro', 'markersize', 10);

        %     prompt = {sprintf('Enter number of frames from %d to delete', k)};
        %     dlg_title = 'Delete future track';
        %     num_lines = 1;
        %     def = {'10'};
        %     answer = inputdlg(prompt,dlg_title,num_lines,def);
            lk=k+eventdata;

            % switch ids
            idx1=(handles.datfile(:,2)==del_id & handles.datfile(:,1)>=k & handles.datfile(:,1) <=lk);
            idx2=(handles.datfile(:,2)==del_id & handles.datfile(:,1) >lk);

            if ~isempty(idx2)
                new_id=max(unique(handles.datfile(:,2)))+1;
                handles.datfile(idx2,2)=new_id;
            end
            handles.datfile(idx1,:)=0;
            instr1=sprintf('[I] Track deleted %d frames from %d', lk-k, k);
            set_instr(instr1, handles);

            guidata(hObject, handles);
        else
            instr1=sprintf('[!]Either the target that you selected is very far from the marked point (Try again)\n');
            instr1=sprintf('%s Or there is no marked target in this frame. Use (a) to add a point.', instr1);
            set_instr(instr1, handles);
        end
    else
        instr1=sprintf('[I] Delete cancelled');
        set_instr(instr1, handles);
    end

end
% --------------------------------------------------------------------
function m_repair_Callback(hObject, eventdata, handles)
% hObject    handle to m_repair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function add_t_Callback(hObject, eventdata, handles)
% hObject    handle to add_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global rp_id
global k
global rp_handle
global ss
global roi_crop
global view_type

if view_type==2 && get(handles.repair_posonly, 'value')
    set(handles.instr, 'string', '[I] Click on the target. (Esc. to cancel)');

    % click on second track (the one you want to switch with the current id)
    [x y button]=ginput(1);
    if button==1        
        plot(x,y, 'k+');
        % if there is no marked positions then copy all from the current
        % set of measurements
        if isempty(handles.frame(k).Mk)
            handles.frame(k).Mk=handles.frame(k).Zk;
        end
        handles.frame(k).Mk=[handles.frame(k).Mk, [x; y; 5]];
        
        
        % proceed if the number of measurements are the same as the number
        % of expected targets
        if ~ss && size(handles.frame(k).Mk,2)==handles.nt 
            ss=1;
            set(handles.startstop, 'string', 'STOP');
            handles=pgn_loopdeloop(handles);
        end
        
        guidata(hObject, handles);
    end
elseif view_type==4

    set(handles.instr, 'string', '[I] Click on the target. (Esc. to cancel)');

    % click on second track (the one you want to switch with the current id)
    [x y button]=ginput(1);
    if button==1
        if exist('rp_handle', 'var')
            if ishandle(rp_handle)
                delete(rp_handle);
            end
        end
        rp_handle=plot(x, y, 'k+');

        id1=find(handles.datfile(:,2)==rp_id & handles.datfile(:,1)<=k);
        [val idx]=max(handles.datfile(id1,1));
        id1=id1(idx);
        X=handles.datfile(id1,:);

        calib=calib2d(roi_crop, handles.sides_cm);
        Zk1=[x; y; 1];

        if ~rp_id
            rp_id=max(unique(handles.datfile(:,2)))+1;
            X(:,2)=rp_id-1;
            X(:,1)=k;
            set_instr(sprintf('[I] No prior fish selected. Creating new id...'), handles);
        end


        if handles.trktype==0 || handles.trktype==1
            % note that the data association for this case is nearest-neighbor
            % since there is only one target and one measurement (marked by the
            % user)
            X=mttkf2d(X, [],  Zk1, 1/handles.fps, calib, 3);
        elseif handles.trktype==2
            X=mttpf2d(X, [],  Zk1, 1/handles.fps, calib);
        end

        curr=X(:,1)==k; % find current time-step
        nzi=find(handles.datfile(:,1)~=0); % find non zero entries
        nzi=nzi(end);
        handles.datfile(nzi+1:nzi+sum(curr),1:size(X,2))=X(curr,:);

        if ~ss
            set_instr(sprintf('[I] Monitor the target fish'), handles);
            ss=1;
            set(handles.startstop, 'string', 'STOP');
            handles=pgn_loopdeloop(handles);
        end

        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function chop_t_Callback(hObject, eventdata, handles)
% hObject    handle to chop_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global k



if view_type==2
    set(handles.instr, 'string', '[I] Click on the target. (Esc. to cancel)');

    % click on second track (the one you want to switch with the current id)
    [x y button]=ginput(1);
    if button==1   
        yp=handles.frame(k).Zk(1,:);
        xp=handles.frame(k).Zk(2,:);
        
        dist=sum(([x; y]*ones(1,numel(xp))-[xp'; yp']).^2);
        [val idx]=min(dist);
        plot(xp(idx), yp(idx), 'ro', 'markersize', 10);
        
        handles.frame(k).Zk(:,idx)=[];
        
        guidata(hObject, handles);
    end
elseif view_type==4

    instr1=sprintf('[I] Click on the track to delete 10 frames from %d frame onwards. (Esc. to cancel)', k);
    set_instr(instr1, handles);


    [x y button]=ginput(1);
    if button==1
        Xk=handles.datfile(handles.datfile(:,1)==k,:);
        pos=Xk(:,3:4);
        [xp yp]=handles.cm2pix(pos(:,1), pos(:,2));
        dist=sum(([x; y]*ones(1,numel(xp))-[xp'; yp']).^2);
        [val idx]=min(dist);

        del_id=Xk(idx, 2);
        plot(xp(idx), yp(idx), 'ro', 'markersize', 10);

        % switch ids
        idx1=(handles.datfile(:,2)==del_id & handles.datfile(:,1)>=k);

        handles.datfile(idx1,:)=0;
        instr1=sprintf('[I] Future instances of this Track deleted');
        set_instr(instr1, handles);

        guidata(hObject, handles);

    else
        instr1=sprintf('[I] Delete cancelled');
        set_instr(instr1, handles);
    end
end

% --- Executes on button press in next_frame.
function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k

k=k+1;
handles=onestep_wrapper(handles); 
set(handles.frmpos, 'value', k);
guidata(hObject, handles);

% --- Executes on button press in prev_frame.
function prev_frame_Callback(hObject, eventdata, handles)
% hObject    handle to prev_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global k

k=k-1;
handles=onestep_wrapper(handles); 
set(handles.frmpos, 'value', k);
guidata(hObject, handles);


% --------------------------------------------------------------------
function help_doc_Callback(hObject, eventdata, handles)
% hObject    handle to help_doc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set_instr('[?] Watch peregrine_help.m4v', handles);
% fprintf('[?] Check your preferences\n');
% fprintf('[?] Did you press save after setting up the thresholds\n');
% fprintf('[?] Did you record the tracks and save?\n');

% --------------------------------------------------------------------
function keypress_short_Callback(hObject, eventdata, handles)
% hObject    handle to keypress_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aboutstr=[];

aboutstr=[aboutstr sprintf('Keyboard shortcuts (the window should be active)... \n')];
aboutstr=[aboutstr sprintf('(Navigation)\n')];
aboutstr=[aboutstr sprintf('Start/Stop [space]\n')];
aboutstr=[aboutstr sprintf('Next/Previous frame [left arrow / right arrow]\n')];
aboutstr=[aboutstr sprintf('\n(Repair mode)\n')];
aboutstr=[aboutstr sprintf('Select target [t]\n')];
aboutstr=[aboutstr sprintf('Switch/update target [s]\n')];
aboutstr=[aboutstr sprintf('Add a point [a]\n')];
aboutstr=[aboutstr sprintf('Delete next few frames [d]\n')];
aboutstr=[aboutstr sprintf('\n(General)\n')];
aboutstr=[aboutstr sprintf('Save [ctrl+s]\n')];


msgbox(aboutstr, 'Keyboard shortcuts');


% --------------------------------------------------------------------
function measure_t_Callback(hObject, eventdata, handles)
% hObject    handle to measure_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rp_id
global k
global roi_crop

set_instr('[I] Click points (at least 3) along the length of the target. Backspace/delete to remove a point. Right click to end)', handles);

% click on second track (the one you want to switch with the current id)
[x y]=getpts;
if numel(x) > 2
    calib=calib2d(roi_crop, handles.sides_cm);
    if rp_id

        id1=(handles.datfile(:,2)==rp_id & handles.datfile(:,1)==k);

        xl=abs(diff(x))*calib.pix2cm(1);
        yl=abs(diff(y))*calib.pix2cm(2);

        fl=sum(sqrt(xl.^2+yl.^2));

        handles.datfile(id1,20)=fl;
        set_instr(sprintf('[I] Measured length=%.1f cm', fl), handles);

        guidata(hObject, handles);
    end
else
    set_instr('[!] Click points (at least 3) along the length of the target. Backspace/delete to remove a point. Right click to end)', handles);
end


% --- Executes on button press in repair_posonly.
function repair_posonly_Callback(hObject, eventdata, handles)
% hObject    handle to repair_posonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of repair_posonly
