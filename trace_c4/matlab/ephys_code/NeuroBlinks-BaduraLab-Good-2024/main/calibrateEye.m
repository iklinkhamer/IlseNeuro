function calibrateEye(obj,event)
%  callback function by video(timer) obj
disp('Delivering puff and saving calibration data.')

cameras = getappdata(0, 'cameras');
metadata = getappdata(0, 'metadata');
src = getselectedsource(cameras(1));

[vid, vid_ts] = getdata(cameras(1), cameras(1).framesAvailable);

% % Set camera to freerun mode so we can preview
% if isprop(src,'FrameStartTriggerSource')
%      src.FrameStartTriggerSource = 'Freerun';
% else
%      src.TriggerSource = 'Freerun';
% end
src.TriggerMode = 'Off';    % Go back to preview mode

% --- save data to root app ---
% Keep data from last trial in memory even if we don't save it to disk
setappdata(0, 'lastvid1', vid);               % Need this for instantReplay
setappdata(0, 'lastmetadata', metadata);
setappdata(0, 'calibration_vid', vid);
setappdata(0, 'calibration_vid_ts', vid_ts);
setappdata(0, 'calibration_metadata', metadata);

fprintf('Video from last trial saved to memory for review.\n')

% metadata.stim.type='None';
                                                 
% --- setting threshold ---H
gui = getappdata(0, 'gui');
gui.eyelidThreshold = ThreshWindowWithPuff;
setappdata(0, 'gui', gui);

drawnow         % Seems necessary to update appdata before returning to calling function

% 
% % Need to allow some time for GUI to draw before we call the lines below
% pause(2)

% Have to do the following 2 lines because we can't call drawhist and
% drawbinary directly from the ThreshWindow opening function since the
% gui struct doesn't exist yet. 
% threshguihandles=guidata(gui.eyelidThreshold);
% ThreshWindowWithPuff('drawbinary',threshguihandles);


