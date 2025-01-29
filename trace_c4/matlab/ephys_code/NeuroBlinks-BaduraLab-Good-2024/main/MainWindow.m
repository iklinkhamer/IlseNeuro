function varargout = MainWindow(varargin)
% MAINWINDOW M-file for MainWindow.fig
%      MAINWINDOW, by itself, creates a new MAINWINDOW or raises the existing
%      singleton*.
%
%      H = MAINWINDOW returns the handle to a new MAINWINDOW or the handle to
%      the existing singleton*.
%
%      MAINWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINWINDOW.M with the given input arguments.
%
%      MAINWINDOW('Property','Value',...) creates a new MAINWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%                    
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainWindow

% Last Modified by GUIDE v2.5 12-Mar-2019 12:21:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainWindow_OpeningFcn, ...
                   'gui_OutputFcn',  @MainWindow_OutputFcn, ...
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


% --- Executes just before MainWindow is made visible.
function MainWindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to MainWindow (see VARARGIN)

% Choose default command line output for MainWindow
handles.output = hObject;

set(handles.axes_eye, 'Color', [240 240 240]/255, 'YAxisLocation', 'right');
set(handles.axes_eye, 'XLim', [-2100 0], 'YLim', [-0.1 1.1])
set(handles.axes_eye, 'XTick', [-3000:500:0], 'box', 'off')
set(handles.axes_eye, 'YTick', [0:0.5:1], 'YTicklabel', {'0' '' '1'})
set(handles.axes_eye, 'FontSize', 8);

% Update handles structure
guidata(hObject,handles)

initializeGUI(hObject, handles)

% --- Executes on button press in pushbutton_StartStopPreview.
function pushbutton_StartStopPreview_Callback(hObject, eventdata, handles)

togglePreview(handles);


function pushbutton_quit_Callback(hObject, eventdata, handles)

gui=getappdata(0,'gui');
config=getappdata(0,'config');

button=questdlg('Are you sure you want to quit?','Quit?');
if ~strcmpi(button,'Yes')
    return
end

ok = stopSession(handles);

% if ~strcmp(config.WHITENOISE_DEVICE_IDS{config.rig}, '')
%     changeWhitenoiseLevel(0);
%     whitenoise_device=getappdata(0,'whitenoise_device');
%     fclose(whitenoise_device);
%     delete(whitenoise_device);
% end

% TODO: Replace these lines with more general loop over all fields of gui?
if isfield(gui, 'one_trial_analysis_gui')
    if ~isempty(gui.one_trial_analysis_gui)
        delete(gui.one_trial_analysis_gui)
    end
end

if isfield(gui, 'camera2gui')
    if ~isempty(gui.camera2gui)
        delete(gui.camera2gui)
    end
end

if isfield(gui, 'maingui')
    if ~isempty(gui.maingui)
        delete(gui.maingui)
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = MainWindow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
varargout{1} = handles.output;


function CamFig_KeyPressFcn(hObject, eventdata, handles)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
switch eventdata.Character
    case '`'
        pushbutton_stim_Callback(hObject, eventdata, handles);
    otherwise
        return
end


function pushbutton_setROI_Callback(hObject, eventdata, handles)

cameras=getappdata(0,'cameras');   
metadata=getappdata(0,'metadata');
gui=getappdata(0,'gui');
config=getappdata(0,'config');

if isfield(metadata.cam(1),'winpos')
    winpos = metadata.cam(1).winpos;
    winpos(1:2) = winpos(1:2) + metadata.cam(1).ROIposition(1:2);
else
    % winpos=[0 0 640 512];
    winpos = config.camera(1).eyelid_roi;
    winpos(1:2) = winpos(1:2) + metadata.cam(1).ROIposition(1:2);
end

% Place rectangle on camera
% h=imrect(handles.cameraAx1,winpos);
h=imellipse(handles.cameraAx1,winpos);

% fcn = makeConstrainToRectFcn('imrect',get(handles.cameraAx1,'XLim'),get(handles.cameraAx1,'YLim'));
fcn = makeConstrainToRectFcn('imellipse',get(handles.cameraAx1,'XLim'),get(handles.cameraAx1,'YLim'));
setPositionConstraintFcn(h,fcn);

% metadata.cam(1).winpos=round(wait(h));
XY=round(wait(h));  % only use for imellipse
metadata.cam(1).winpos=round(getPosition(h));
metadata.cam(1).winpos(1:2)=metadata.cam(1).winpos(1:2)-metadata.cam(1).ROIposition(1:2);
metadata.cam(1).mask=createMask(h, handles.pwin(1));

wholeframe=getsnapshot(cameras(1));
binframe=im2bw(wholeframe,metadata.cam(1).thresh);
eyeframe=binframe.*metadata.cam(1).mask;
metadata.cam(1).pixelpeak=sum(sum(eyeframe));

hp=findobj(handles.cameraAx1,'Tag','roipatch');
delete(hp)
% handles.roipatch=patch([xmin,xmin+width,xmin+width,xmin],[ymin,ymin,ymin+height,ymin+height],'g','FaceColor','none','EdgeColor','g','Tag','roipatch');
% XY=getVertices(h);
delete(h);
handles.roipatch=patch(XY(:,1),XY(:,2),'g','FaceColor','none','EdgeColor','g','Tag','roipatch');
handles.XY=XY;

setappdata(0,'metadata',metadata);
guidata(gui.maingui,handles)

drawnow         % Seems necessary to update appdata before returning to calling function


function pushbutton_calibrateEye_Callback(hObject, eventdata, handles)

microControllerVariablesEnum;

metadata=getappdata(0, 'metadata'); 
metadata.cam(1).cal = 1;

metadata = refreshParams(handles, metadata);
ok = uploadParams(metadata);

if ~ok
    msgbox('Problem communicating with Microcontroller, aborting calibration')
    return
end

cameras=getappdata(0,'cameras');
cameras(1).TriggerRepeat = 0;
cameras(1).StopFcn=@calibrateEye;   % this will be executed after timer stop
flushdata(cameras(1));         % Remove any data from buffer before triggering

% Set camera to hardware trigger mode
src = getselectedsource(cameras(1));
% src.FrameStartTriggerSource = 'Line1';
src.TriggerSource='Line1'; %turned on by NHE. IK note 2-11-23 Very important that this line remains.
% src.TriggerMode = 'On';     % Normally set it to hardware trigger turned % off by NHE. IK note 2-11-23 Very important that this line remains off.
% Placed here by MvdDoe, supresses warning that cloggs up the command window
warning('off','imaq:gentl:hardwareTriggerTriggerModeOff')
% TODO: REMOVE THE FOLLOWING LINE BEFORE PRODUCTION
%tic
start(cameras(1));
% src.TriggerMode = 'Off';    % For now override hardware trigger so we don't need Arduino connected
metadata.cam(1).cal = 0;

metadata.ts(2) = etime(clock, datevec(metadata.ts(1)));
%toc

% --- trigger via microController --
microController=getappdata(0, 'microController');
fwrite(microController, uController.RUN, 'uint8');

setappdata(0, 'metadata', metadata);

drawnow         % Seems necessary to update appdata before returning to calling function

% --- Executes on button press in togglebutton_tgframerate.
function togglebutton_tgframerate_Callback(hObject, eventdata, handles)

    % Not currently implemented -- need to rewrite

% cameras=getappdata(0,'cameras');
% metadata=getappdata(0,'metadata');

% src = getselectedsource(cameras(1));

% if get(hObject,'Value')
%     % Turn on high frame rate mode
%     metadata.cam(1).camera_ROIposition=max(metadata.cam(1).winpos+[-10 0 20 0],[0 0 0 0]);
%     cameras(1).ROIposition=metadata.cam(1).camera_ROIposition;
% %     metadata.cam(1).frame_rate=500;
%     src.ExposureTimeAbs = 1900;
% %     src.AllGainRaw=metadata.cam(1).init_AllGainRaw+round(20*log10(metadata.cam(1).init_ExposureTime/src.ExposureTimeAbs));
%     % --- size fit for roi and mask ----
%     vidroi_x=metadata.cam(1).camera_ROIposition(1)+[1:metadata.cam(1).camera_ROIposition(3)];
%     vidroi_y=metadata.cam(1).camera_ROIposition(2)+[1:metadata.cam(1).camera_ROIposition(4)];
%     metadata.cam(1).mask = metadata.cam(1).mask(vidroi_y, vidroi_x);
%     metadata.cam(1).winpos(1:2)=metadata.cam(1).winpos(1:2)-metadata.cam(1).camera_ROIposition(1:2);
% else
%     % Turn off high frame rate mode
%     camera.ROIposition=metadata.cam(1).fullsize;
% %     metadata.cam(1).frame_rate=200;
%     src.ExposureTimeAbs = metadata.cam(1).init_ExposureTime;
% %     src.AllGainRaw=metadata.cam(1).init_AllGainRaw;
%     % --- size fit for roi and mask ----
%     mask0=metadata.cam(1).mask; 
%     s_mask0=size(mask0);
%     metadata.cam(1).mask = false(metadata.cam(1).fullsize([4 3]));
%     metadata.cam(1).mask(metadata.cam(1).camera_ROIposition(2)+[1:s_mask0(1)], metadata.cam(1).camera_ROIposition(1)+[1:s_mask0(2)])=mask0;
%     metadata.cam(1).winpos(1:2)=metadata.cam(1).winpos(1:2)+metadata.cam(1).camera_ROIposition(1:2);
%     metadata.cam(1).camera_ROIposition=metadata.cam(1).fullsize;
% end

% % pushbutton_StartStopPreview_Callback(handles.pushbutton_StartStopPreview, [], handles)
% % pause(0.02)
% % pushbutton_StartStopPreview_Callback(handles.pushbutton_StartStopPreview, [], handles)

% setappdata(0,'camera',camera);
% setappdata(0,'metadata',metadata);

drawnow         % Seems necessary to update appdata before returning to calling function

function checkbox_record_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.checkbox_record,'BackgroundColor',[0 1 0]); % green
else
    set(handles.checkbox_record,'BackgroundColor',[1 0 0]); % red
end


function pushbutton_instantreplay_Callback(hObject, eventdata, handles)
instantReplay(getappdata(0,'lastvid1'),getappdata(0,'lastmetadata'));


function toggle_continuous_Callback(hObject, eventdata, handles)
if get(hObject,'Value'),
    checkContext(handles) %turned off NHE
    set(hObject,'String','Pause Continuous')
    set(handles.trialtimecounter,'Visible','On')
else
    set(hObject,'String','Start Continuous')
    set(handles.trialtimecounter,'Visible','Off')
end


function pushbutton_singleTrial_Callback(hObject, eventdata, handles)
checkContext(handles)
startTrial(handles)

function popupmenu_stimtype_Callback(hObject, eventdata, handles)
% --- updating metadata ---
metadata=getappdata(0,'metadata');
val=get(hObject,'Value');
str=get(hObject,'String');
metadata.stim.type=str{val};
setappdata(0,'metadata',metadata);


% ------ highlight for uipanel -----
set(handles.uipanel_puff,'BackgroundColor',[240 240 240]/255);
set(handles.uipanel_conditioning,'BackgroundColor',[240 240 240]/255);
switch lower(metadata.stim.type)
    case 'puff'
        set(handles.uipanel_puff,'BackgroundColor',[225 237 248]/255); % light blue
    case 'conditioning'
        set(handles.uipanel_conditioning,'BackgroundColor',[225 237 248]/255); % light blue
end   

drawnow         % Seems necessary to update appdata before returning to calling function

function togglebutton_stream_Callback(hObject, eventdata, handles)

if get(hObject,'Value'),
    startStreaming(handles)
else
    stopStreaming(handles)
end


function pushbutton_params_Callback(hObject, eventdata, handles)
ParamsWindow


function pushbutton_eyelidTraceViewer_Callback(hObject, eventdata, handles)
gui=getappdata(0,'gui');
config=getappdata(0,'config');
gui.one_trial_analysis_gui = OneTrialAnalysisWindow;

set(gui.one_trial_analysis_gui, 'units', 'pixels')
current_pos = get(gui.one_trial_analysis_gui, 'position');
set(gui.one_trial_analysis_gui, 'position', [config.pos_analysiswindow{config.rig}, current_pos(3:4)]); % mod 2 because we have 2 rigs per computer/monitor

setappdata(0,'gui',gui);

function uipanel_SessionMode_SelectionChangeFcn(hObject, eventdata, handles)

metadata = getappdata(0, 'metadata');
cameras = getappdata(0, 'cameras');

switch get(eventdata.NewValue, 'Tag') % Get Tag of selected object.
    case 'togglebutton_NewSession'
        dlgans = inputdlg({'Enter session number'}, 'Create');
        if isempty(dlgans) 
            ok=0;
        elseif isempty(dlgans{1})
            ok=0;
        else
            ok=1;       
            metadata.session=str2num(dlgans{1});
        end

        setappdata(0, 'metadata', metadata)

        if ok
%             try
                % Consider passing 'config' struct instead of handles and storing session number in 'config'
                ok = startSession(handles);
%             catch ME
%                 set(eventdata.NewValue,'Value',0);
%                 set(eventdata.OldValue,'Value',1);
%                 set(handles.uipanel_SessionMode,'SelectedObject',eventdata.OldValue);
%                 throw(ME)
%             end
        end

    case 'togglebutton_StopSession'
        button=questdlg('Are you sure you want to stop this session?', 'Stop session?', ...
        'Yes and compress videos', 'Yes and DON''T compress videos', 'No', 'Yes and compress videos');
        
        switch button
            
            case 'Yes and compress videos'
                metadata.session=0;     
                ok=1;
                stopSession(handles);
                
                makeCompressedVideos(metadata.folder,1);
                
            case 'Yes and DON''T compress videos'
                metadata.session=0;     
                ok=1;
                stopSession(handles);
                
            otherwise
                ok=0;

            setappdata(0, 'metadata', metadata)
                
        end
    otherwise
        warndlg('There is something wrong with the mode selection callback','Mode Select Problem!')
        return
end

if ok
    set(eventdata.NewValue,'Value',1);
    set(eventdata.OldValue,'Value',0);
    set(handles.uipanel_SessionMode,'SelectedObject',eventdata.NewValue);
else
    set(eventdata.NewValue,'Value',0);
    set(eventdata.OldValue,'Value',1);
    set(handles.uipanel_SessionMode,'SelectedObject',eventdata.OldValue);
    return
end

set(handles.text_SessionName, 'String', metadata.session);

drawnow         % Seems necessary to update appdata before returning to calling function


function pushbutton_opentable_Callback(hObject, eventdata, handles)
% paramtable.data=get(handles.uitable_params,'Data');
% paramtable.randomize=get(handles.checkbox_random,'Value');
% % paramtable.tonefreq=str2num(get(handles.edit_tone,'String'));
% % if length(paramtable.tonefreq)<2, paramtable.tonefreq(2)=0; end
% setappdata(0,'paramtable', paramtable);

drawnow         % Seems necessary to update appdata before returning to calling function

gui=getappdata(0,'gui');
trialtablegui=TrialTable;
% movegui(trialtablegui,[gui.pos_mainwin(1)+gui.size_mainwin(1)+20 gui.pos_mainwin(2)])



% --- Executes on button press in checkbox_random.
function checkbox_random_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_random


function edit_tone_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_tone as text
%        str2double(get(hObject,'String')) returns contents of edit_tone as a double


% --- Executes during object creation, after setting all properties.
function edit_tone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function popupmenu_stimtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_stimtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pretime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pretime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_pretime as text
%        str2double(get(hObject,'String')) returns contents of edit_pretime as a double


% --- Executes during object creation, after setting all properties.
function edit_pretime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pretime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in pushbutton7.
% function pushbutton7_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton7 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)


function edit_ITI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ITI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_ITI as text
%        str2double(get(hObject,'String')) returns contents of edit_ITI as a double


% --- Executes during object creation, after setting all properties.
function edit_ITI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ITI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_puffdur_Callback(hObject, eventdata, handles)
% hObject    handle to edit_puffdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_puffdur as text
%        str2double(get(hObject,'String')) returns contents of edit_puffdur as a double

value = str2double(hObject.String);
if value > 50
    warndlg(sprintf('WARNING: You''re setting a higher puff duration (%d ms) than is typical! Consider changing', value))
end
   


% --- Executes during object creation, after setting all properties.
function edit_puffdur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_puffdur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_SessionName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_SessionName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





function edit_stabletime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stabletime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stabletime as text
%        str2double(get(hObject,'String')) returns contents of edit_stabletime as a double


% --- Executes during object creation, after setting all properties.
function edit_stabletime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stabletime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stableeye_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stableeye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stableeye as text
%        str2double(get(hObject,'String')) returns contents of edit_stableeye as a double


% --- Executes during object creation, after setting all properties.
function edit_stableeye_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stableeye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_eyelidthresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eyelidthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eyelidthresh as text
%        str2double(get(hObject,'String')) returns contents of edit_eyelidthresh as a double


% --- Executes during object creation, after setting all properties.
function edit_eyelidthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eyelidthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_verbose.
function checkbox_verbose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

config = getappdata(0, 'config');

% Hint: get(hObject,'Value') returns toggle state of checkbox_verbose
if hObject.Value == 0
    config.verbose = 0;
    % Suppress some warnings that are annoying
    warning('off', 'imaq:gentl:immedOrManTriggerTurningTriggerModeOff')
else
    config.verbose = 1;
    warning('on', 'imaq:gentl:immedOrManTriggerTurningTriggerModeOff')
end

setappdata(0, 'config', config);

function edit_posttime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_posttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_posttime as text
%        str2double(get(hObject,'String')) returns contents of edit_posttime as a double


% --- Executes during object creation, after setting all properties.
function edit_posttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_posttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_abort.
function pushbutton_abort_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If camera gets hung up for any reason, this button can be pressed to
% reset it.

abortCameraAcquisition;




function edit_ITI_rand_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ITI_rand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ITI_rand as text
%        str2double(get(hObject,'String')) returns contents of edit_ITI_rand as a double


% --- Executes during object creation, after setting all properties.
function edit_ITI_rand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ITI_rand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_StopAfterTrial_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StopAfterTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StopAfterTrial as text
%        str2double(get(hObject,'String')) returns contents of edit_StopAfterTrial as a double


% --- Executes during object creation, after setting all properties.
function edit_StopAfterTrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StopAfterTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_loadParams.
function pushbutton_loadParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

checkContext(handles)
loadParamTable(handles);



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function checkbox_verbose_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

config = getappdata(0, 'config');

% Hint: get(hObject,'Value') returns toggle state of checkbox_verbose
if config.verbose == 0
    hObject.Value = 0;   
else
    hObject.Value = 1;
end


% --- Executes during object creation, after setting all properties.
function checkbox_random_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
config = getappdata(0, 'config');

% Hint: get(hObject,'Value') returns toggle state of checkbox_verbose
if config.paramtable.randomize == 0
    hObject.Value = 0;
    
else
    hObject.Value = 1;
end


% --- Executes when user attempts to close CamFig.
function CamFig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CamFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
msgbox('This window cannot be closed during an experiment')


% --- Executes when entered data in editable cell(s) in uitable_params.
function uitable_params_CellEditCallback(hObject, eventdata, handles)
    % hObject    handle to uitable_params (see GCBO)
    % eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
    %	Indices: row and column indices of the cell(s) edited
    %	PreviousData: previous data for the cell(s) edited
    %	EditData: string(s) entered by the user
    %	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
    %	Error: error string when failed to convert EditData to appropriate value for Data
    % handles    structure with handles and user data (see GUIDATA)
    
    config = getappdata(0, 'config');

    config.paramtable.data = get(handles.uitable_params, 'Data');
    
    setappdata(0, 'config', config)



% --- Executes on selection change in popupmenu_whitenoise_control.
function popupmenu_whitenoise_control_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_whitenoise_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_whitenoise_control contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_whitenoise_control

selection = hObject.String{hObject.Value};

switch selection
    case 'Off'
        changeWhitenoiseLevel(0);
    case 'On'
        changeWhitenoiseLevel(1);
    case 'Trial table'
        changeWhitenoiseLevel(0);
end


% --- Executes during object creation, after setting all properties.
function popupmenu_whitenoise_control_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_whitenoise_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
