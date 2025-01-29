function varargout = Camera2Window(varargin)
% CAMERAWINDOW MATLAB code for cameraWindow.fig
%      CAMERAWINDOW, by itself, creates a new CAMERAWINDOW or raises the existing
%      singleton*.
%
%      H = CAMERAWINDOW returns the handle to a new CAMERAWINDOW or the handle to
%      the existing singleton*.
%
%      CAMERAWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAMERAWINDOW.M with the given input arguments.
%
%      CAMERAWINDOW('Property','Value',...) creates a new CAMERAWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cameraWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cameraWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cameraWindow

% Last Modified by GUIDE v2.5 25-Oct-2018 14:25:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Camera2Window_OpeningFcn, ...
                   'gui_OutputFcn',  @Camera2Window_OutputFcn, ...
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


% --- Executes just before camera2Window is made visible.
function Camera2Window_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to camera2Window (see VARARGIN)


% Choose default command line output for camera2Window
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes camera2Window wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Camera2Window_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_instantReplay.
function pushbutton_instantReplay_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_instantReplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cam2 = getappdata(0,'cam2');

instantReplay(getappdata(0,'lastvid2'),getappdata(0,'lastmetadata'))



% % --- Executes on button press in togglebutton_startStopPreview.
% function togglebutton_startStopPreview_Callback(hObject, eventdata, handles)
% % hObject    handle to togglebutton_startStopPreview (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% cameras = getappdata(0,'cameras');

% if length(cameras) < 2
%     return
% end

% if hObject.Value == 1
%     set(hObject,'String','Stop Preview')
%     pwin=image(zeros(480,640), 'Parent',handles.cameraAx2);
%     preview(cameras(2),pwin)
% else
%     set(hObject,'String','Start Preview')
%     stoppreview(cameras(2))
% end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% config = getappdata(0, 'config');
% 
% % Hint: delete(hObject) closes the figure
% if config.running
%     msgbox('This window cannot be closed during an experiment')
% else
%     delete(hObject);
% end

msgbox('This window cannot be closed during an experiment')
