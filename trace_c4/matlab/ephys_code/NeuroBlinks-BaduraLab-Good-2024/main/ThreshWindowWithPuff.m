function varargout = ThreshWindowWithPuff(varargin)
% THRESHWINDOWWITHPUFF MATLAB code for ThreshWindowWithPuff.fig
%      THRESHWINDOWWITHPUFF, by itself, creates a new THRESHWINDOWWITHPUFF or raises the existing
%      singleton*.
%
%      H = THRESHWINDOWWITHPUFF returns the handle to a new THRESHWINDOWWITHPUFF or the handle to
%      the existing singleton*.
%
%      THRESHWINDOWWITHPUFF('CALLBACK', hObject, eventData, handles, ...) calls the local
%      function named CALLBACK in THRESHWINDOWWITHPUFF.M with the given input arguments.
%
%      THRESHWINDOWWITHPUFF('Property', 'Value', ...) creates a new THRESHWINDOWWITHPUFF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThreshWindowWithPuff_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThreshWindowWithPuff_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThreshWindowWithPuff

% Last Modified by GUIDE v2.5 09-Dec-2013 16:27:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThreshWindowWithPuff_OpeningFcn, ...
                   'gui_OutputFcn',  @ThreshWindowWithPuff_OutputFcn, ...
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


% --- Executes just before ThreshWindowWithPuff is made visible.
function ThreshWindowWithPuff_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThreshWindowWithPuff (see VARARGIN)

metadata = getappdata(0, 'calibration_metadata');

% Choose default command line output for ThreshWindowWithPuff
handles.output = hObject;

handles.x1 = floor(metadata.cam(1).winpos(1)) + 1;
handles.x2 = floor(metadata.cam(1).winpos(1) + metadata.cam(1).winpos(3));
handles.y1 = floor(metadata.cam(1).winpos(2)) + 1;
handles.y2 = floor(metadata.cam(1).winpos(2) + metadata.cam(1).winpos(4));

set(handles.edit_eyelidThresh, 'String', num2str(round(metadata.cam(1).thresh * 255)));
set(handles.slider_eyelidThresh, 'Value', round(metadata.cam(1).thresh * 255));

% Update handles structure
guidata(hObject, handles);

drawbinary(handles, hObject);

% UIWAIT makes ThreshWindowWithPuff wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% edit_eyelidThresh_Callback(handles.edit_eyelidThresh, 0, handles)


% --- Outputs from this function are returned to the command line.
function varargout = ThreshWindowWithPuff_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_eyelidThresh_Callback(hObject, eventdata, handles)
% hObject    handle to slider_eyelidThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject, 'Value') returns position of slider
%        get(hObject, 'Min') and get(hObject, 'Max') to determine range of slider

metadata = getappdata(0, 'calibration_metadata');

% Alternatively can use "thresh = graythresh(image)" to get threshold with
% Otsu's method
threshold_slider = round(get(hObject, 'Value')) / 255;
if threshold_slider > 0 && threshold_slider < 1
    metadata.cam(1).thresh = round(get(hObject, 'Value')) / 255;
    set(handles.edit_eyelidThresh, 'String', num2str(round(get(hObject, 'Value'))));
else
    disp('Wrong number was returned by slider...')
    set(handles.slider_eyelidThresh, 'Value', round(metadata.cam(1).thresh * 255));
end
% Have to stop auto refresh timer before drawing or else we get an error
% when the timer has a callback
% stop(handles.timer)
setappdata(0, 'calibration_metadata', metadata);
drawbinary(handles);


% --- Executes during object creation, after setting all properties.
function slider_eyelidThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_eyelidThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', [.9 .9 .9]);
end



function edit_eyelidThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eyelidThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject, 'String') returns contents of edit_eyelidThresh as text
%        str2double(get(hObject, 'String')) returns contents of edit_eyelidThresh as a double

metadata = getappdata(0, 'calibration_metadata');

% Alternatively can use "thresh = graythresh(image)" to get threshold with
% Otsu's method
metadata.cam(1).thresh = str2double(get(handles.edit_eyelidThresh, 'String'))/255;
set(handles.slider_eyelidThresh, 'Value', round(metadata.cam(1).thresh * 255));

setappdata(0, 'calibration_metadata', metadata);

drawbinary(handles);


% --- Executes during object creation, after setting all properties.
function edit_eyelidThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eyelidThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
% calculateEyelidCalibration;      % Don't do calibration if user closes window without pressing enter

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);
calculateEyelidCalibration;


% --- Executes on button press in pushbutton_only_thresh.
function pushbutton_only_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_only_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

metadata = getappdata(0, 'calibration_metadata');
fprintf('thresh = %d.\n', round(metadata.cam(1).thresh * 255))
delete(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% user defined functions %%%%%%%%%%

function drawbinary(handles, varargin)

% Load objects from root app data
gui = getappdata(0, 'gui');  % Load global handles
metadata = getappdata(0, 'calibration_metadata');
vid = getappdata(0, 'calibration_vid');

% This is a hack to deal with the case where drawbinary has to be called
% before figure has been fully created (in Opening function)
if length(varargin) > 0 
    hf = varargin{1};
else
    hf = gui.eyelidThreshold;
end

% --- eyelid trace --
[trace, t] = vid2eyetrace(vid, metadata, metadata.cam(1).thresh);
ind_t = find(t < 0.2);
[y_max, ind_max1] = max(trace(ind_t));  [y_min, ind_min1] = min(trace(ind_t));
ind_max = ind_t(ind_max1);  ind_min = ind_t(ind_min1);  
set(0, 'CurrentFigure', hf)

% --- for axex_binary: eye opened ---
wholeframe = vid(:, :, 1, ind_min(1));
% roi = wholeframe(handles.y1:handles.y2, handles.x1:handles.x2);
% Had to revise this to work with elliptical ROI
roi = wholeframe .* uint8(metadata.cam(1).mask);
binframe = im2bw(roi(handles.y1:handles.y2, handles.x1:handles.x2), metadata.cam(1).thresh);
handles.binimage = imshow(binframe, 'Parent', handles.axes_binary);

set(hf, 'CurrentAxes', handles.axes_hist)
imhist(roi(metadata.cam(1).mask))

% --- for axex_binary: eye closeed ---
wholeframe = vid(:, :, 1, ind_max(1));
% roi = wholeframe(handles.y1:handles.y2, handles.x1:handles.x2);
% Had to revise this to work with elliptical ROI
roi = wholeframe .* uint8(metadata.cam(1).mask);
binframe = im2bw(roi(handles.y1:handles.y2, handles.x1:handles.x2), metadata.cam(1).thresh);
handles.binimage2 = imshow(binframe, 'Parent', handles.axes_binary2);

set(hf, 'CurrentAxes', handles.axes_hist2)
imhist(roi(metadata.cam(1).mask))

% --- for axex_trace: eyelid trace ---
% set(hf, 'CurrentAxes', handles.axes_trace)
plot(handles.axes_trace, t, trace)
set(handles.axes_trace, 'xlim', [t(1) t(end)])
