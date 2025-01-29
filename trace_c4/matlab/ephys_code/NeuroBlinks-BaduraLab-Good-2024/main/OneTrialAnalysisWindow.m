function varargout = OneTrialAnalysisWindow(varargin)
% OneTrialAnalysisWindow MATLAB code for OneTrialAnalysisWindow.fig
%      OneTrialAnalysisWindow, by itself, creates a new OneTrialAnalysisWindow or raises the existing
%      singleton*.
%
%      H = OneTrialAnalysisWindow returns the handle to a new OneTrialAnalysisWindow or the handle to
%      the existing singleton*.
%
%      OneTrialAnalysisWindow('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OneTrialAnalysisWindow.M with the given input arguments.
%
%      OneTrialAnalysisWindow('Property','Value',...) creates a new OneTrialAnalysisWindow or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OneTrialAnalysisWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OneTrialAnalysisWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OneTrialAnalysisWindow

% Last Modified by GUIDE v2.5 25-Mar-2013 15:19:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @OneTrialAnalysisWindow_OpeningFcn, ...
    'gui_OutputFcn',  @OneTrialAnalysisWindow_OutputFcn, ...
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


% --- Executes just before OneTrialAnalysisWindow is made visible.
function OneTrialAnalysisWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OneTrialAnalysisWindow (see VARARGIN)

% Choose default command line output for OneTrialAnalysisWindow
handles.output = hObject;
handles.tnum_prev=NaN;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OneTrialAnalysisWindow wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OneTrialAnalysisWindow_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
t=timerfind('Name','analysisUpdateTimer');
if ~isempty(t)
    stop(t);
end
delete(t)
catch
    disp('Error in closing fig.')
end
% Hint: delete(hObject) closes the figure
delete(hObject);





% --- Executes on slider movement.
function slider_trialnum_Callback(hObject, eventdata, handles)
% hObject    handle to slider_trialnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if strcmp(get(handles.edit_trialnum,'String'),'Trial Num')
    set(handles.edit_trialnum,'String','0');
end
t_num=round(str2double(get(handles.edit_trialnum,'String'))+get(hObject,'Value'));

trials=getappdata(0,'trials');
if isfield(trials,'eye')
    tnum_max=length(trials.eye);
else
    tnum_max=0;
end

% set(handles.text_trialnum,'String',sprintf('Trial Viewer (1-%d)',tnum_max));

if t_num>tnum_max
    t_num=tnum_max;
elseif t_num<1
    t_num=1;
end
set(handles.edit_trialnum,'String',num2str(t_num)),
set(hObject,'Value',0),

if t_num~=handles.tnum_prev | t_num==tnum_max,
    drawOneEyelid(handles,t_num);   
end
handles.tnum_prev=t_num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_trialnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_trialnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',0);


function edit_trialnum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trialnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trialnum as text
%        str2double(get(hObject,'String')) returns contents of edit_trialnum as a double

t_num=round(str2double(get(handles.edit_trialnum,'String')));
trials=getappdata(0,'trials');
if isfield(trials,'eye')
    tnum_max=length(trials.eye);
else
    tnum_max=0;
end
if t_num>tnum_max || t_num<1
    if t_num<1
        t_num=1;  
    else
        t_num=tnum_max;
    end
    set(handles.edit_trialnum,'String',num2str(t_num)),
end

if t_num~=handles.tnum_prev | t_num==tnum_max,
    drawOneEyelid(handles,t_num);   
end
handles.tnum_prev=t_num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_trialnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trialnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% user difined functions %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawOneEyelid(handles,t_num)

trials=getappdata(0,'trials');
if ~isfield(trials,'eye'), return,  end
tnum_max=length(trials.eye);
str_tl={[sprintf('Trial Viewer (1-%d)',tnum_max)]};

% ------- for eye -----
cla(handles.axes_eyelid)

set(handles.axes_eyelid,'ylim',[-0.15 1.20], 'ytick',[0:0.5:1], 'box', 'off','tickdir','out', 'NextPlot', 'add')

plot(handles.axes_eyelid, [-1 1]*1000, [0 0],'k:')
plot(handles.axes_eyelid, [-1 1]*1000, [1 1],'k:')

xlim1=[trials.eye(t_num).time(1) trials.eye(t_num).time(end)];
plotOneEyelid(handles.axes_eyelid, t_num);

% text(xlim1*[0.33;0.67], -0.46, str_tl)
set(handles.axes_eyelid,'xlim',xlim1,'xtick',[-400:200:1000])
set(handles.axes_eyelid,'color',[240 240 240]/255);

