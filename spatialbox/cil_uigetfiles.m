function varargout = cil_uigetfiles(varargin)
% CIL_UIGETFILES M-file for cil_uigetfiles.fig
%      CIL_UIGETFILES, by itself, creates a new CIL_UIGETFILES or raises the existing
%      singleton*.
%
%      H = CIL_UIGETFILES returns the handle to a new CIL_UIGETFILES or the handle to
%      the existing singleton*.
%
%      CIL_UIGETFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CIL_UIGETFILES.M with the given input arguments.
%
%      CIL_UIGETFILES('Property','Value',...) creates a new CIL_UIGETFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cil_uigetfiles_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cil_uigetfiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)"
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Read: http://www.mathworks.com/support/solutions/data/33916.html
%       http://www.mathworks.com/support/solutions/data/29133.html

% Edit the above text to modify the response to help cil_uigetfiles
% Last Modified by Alex Liberzon v2.5 18-Apr-2004 20:30:03
% CIL Flow llc. (c) Copyright 2004

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cil_uigetfiles_OpeningFcn, ...
    'gui_OutputFcn',  @cil_uigetfiles_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    %     str2func(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cil_uigetfiles is made visible.
function cil_uigetfiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cil_uigetfiles (see VARARGIN)

% Choose default command line output for cil_uigetfiles
handles = guihandles;
handles.output = hObject;
handles.path = cd;
handles.state3d = 0;

% Update handles structure
guidata(hObject, handles);
    update_gui(hObject,[],handles);
% % 
% % % UIWAIT makes cil_uigetfiles wait for user response (see %uiresume)
uiwait(handles.fig);


% --- Outputs from this function are returned to the command line.
function varargout = cil_uigetfiles_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% if nargout
if isfield(handles,'filenames')
    varargout{1} = handles.filenames;
    varargout{2} = handles.path;
    varargout{3} = handles.dT*handles.step;
    varargout{4} = handles.scale;
    varargout{5} = handles.state3d;
    
else
    varargout{1} = {};
    varargout{2} = {};
    varargout{3} = NaN;
    varargout{4} = NaN;
    varargout{5} = NaN;
end
close(handles.fig);


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Returns the names of the selected files
handles.list_entries = get(handles.listbox_files,'String');
handles.index_selected = get(handles.listbox_files,'Value');
if isempty(handles.index_selected) 
    errordlg('You must select file before pressing Load','Incorrect Selection','modal')
end

% if ~isfield(handles,'step')
handles.step = str2double(get(handles.edit_step,'String'));
handles.dT = str2double(get(handles.edit_dT,'String'));
handles.scale = str2double(get(handles.edit_scale,'String'));
% end

switch length(handles.index_selected)
    case {1}                        % only the first file is selected,
        % pick all files from it up to last
        % How many files are selected
        
        index = handles.index_selected:handles.step:length(handles.list_entries);
        % [handles.filenames{1:length(index),1}] = deal(handles.list_entries(index));
        handles.filenames = handles.list_entries(index);
    case {2}
        % two files are selected, first and last
        index = min(handles.index_selected):handles.step:max(handles.index_selected);
        %         [handles.filenames{1:length(index),1}] = deal(handles.list_entries(index));
        handles.filenames = handles.list_entries(index);
        
    otherwise 
        % many files are selected, no step is involved, just read
        % all of them
        % [handles.filenames{1:length(handles.index_selected),1}] = deal(handles.list_entries(handles.index_selected));
        handles.filenames = handles.list_entries(handles.index_selected);
        
end
guidata(hObject,handles);
uiresume(handles.fig);




% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.fig);


% --- Executes during object creation, after setting all properties.
function listbox_files_CreateFcn(hObject, eventdata, handles)

% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end
% if ~isfield(handles,'path')
%     handles.path = cd;
% end
% 
% if get(handles.check3d,'Value') == 1
%     handles.files = dir(fullfile(handles.path,'*.v3d')); %     tmp2
% else
%     handles.files = dir(fullfile(handles.path,'*.vec'));
% end;  
% 
% list = dir(handles.path);
% ind = find(cat(1,list.isdir));
% set(handles.listbox_files,'String',{list(ind).name,handles.files.name});
% guidata(hObject,handles)

% --- Executes on selection change in listbox_files.
function listbox_files_Callback(hObject, eventdata, handles)

if strcmp(get(handles.fig,'SelectionType'),'open') % If double click
    index_selected = get(handles.listbox_files,'Value');
    file_list = get(handles.listbox_files,'String');
    filename = file_list{index_selected};       % Item selected in list box
%     keyboard
    if  isdir([handles.path,filesep,filename]) % If directory
        if index_selected == 2
            handles.path = handles.path(1:max(findstr(handles.path,filesep)-1));
        elseif index_selected > 2
            handles.path = [handles.path,filesep,filename];
        end
        guidata(handles.fig,handles);
        update_gui(handles.fig,[],handles);
    end
end




% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
handles = guihandles;
if ~isfield(handles,'path')
    handles.path = cd;
end
set(hObject,'String',handles.path);
guidata(hObject,handles);



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double

tmp = get(hObject,'String');
if isdir(tmp)
    handles.path = tmp; 
else
    set(hObject,'String',handles.path);
end
guidata(hObject,handles);
update_gui(hObject, [], handles);

% --- Executes on button press in pushbutton_cd.
function pushbutton_cd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = uigetdir(handles.path);
if tmp > 0 & isdir(tmp), handles.path = tmp; end
set(handles.edit_path,'String',handles.path);
guidata(hObject, handles);
update_gui(hObject, [], handles);


function edit_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step as text
%        str2double(get(hObject,'String')) returns contents of edit_step as a double
tmp = str2double(get(hObject,'String'));
if isnan(tmp) | floor(tmp) ~= tmp
    set(hObject,'String','1');
else
    handles.step = tmp;
end
guidata(hObject,handles);

function edit_dT_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dT as text
%        str2double(get(hObject,'String')) returns contents of edit_dT as a double
if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','1');
end
handles.dT = str2double(get(hObject,'String'));
guidata(hObject,handles);


function edit_scale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scale as text
%        str2double(get(hObject,'String')) returns contents of edit_scale as a double
if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','1');
end
handles.scale = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles = guihandles(hObject);
movegui(hObject,'northwest')


% --- Executes on button press in pushbutton_updir.
function pushbutton_updir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_updir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = handles.path;
s = s(1:max(3,max(findstr(s,filesep))-1));
if exist(s,'dir')
    handles.path = s;
end
guidata(hObject,handles);
update_gui(hObject,[],handles);

function update_gui(hObject, eventdata, handles)
% Self made UPDATE GUI function
% %keyboard
if isdir(handles.path)
    set(handles.edit_path,'String',handles.path);
    if get(handles.check3d,'Value') == 1
        handles.files = dir(fullfile(handles.path,'*.v3d')); %     tmp2
    else
        handles.files = dir(fullfile(handles.path,'*.vec'));
    end;   
else
    handles.path = cd;
    set(handles.edit_path,'String',handles.path);
    if get(handles.check3d,'Value') == 1
        handles.files = dir(fullfile(handles.path,'*.v3d')); %     tmp2
    else
        handles.files = dir(fullfile(handles.path,'*.vec'));
    end;   
end
list = dir(handles.path);
ind = find(cat(1,list.isdir));
set(handles.fig,'SelectionType','normal');
set(handles.listbox_files,'String',{list(ind).name,handles.files.name},'Value',1);    
guidata(handles.fig, handles);



% --- Executes during object creation, after setting all properties.
function edit_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
tmp = str2double(get(hObject,'String'));
if isnan(tmp) | floor(tmp) ~= tmp
    set(hObject,'String','1');
else
    handles.step = tmp;
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_dT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','1');
end
handles.dT = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


if isnan(str2double(get(hObject,'String')))
    set(hObject,'String','1');
end
handles.scale = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes on button press in check3d.
function check3d_Callback(hObject, eventdata, handles)
% hObject    handle to check3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check3d
handles.state3d  = 1;
guidata(hObject,handles);
update_gui(hObject, [], handles);

