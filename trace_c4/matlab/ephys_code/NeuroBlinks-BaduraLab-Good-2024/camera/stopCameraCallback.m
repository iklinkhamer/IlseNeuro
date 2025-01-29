function stopCameraCallback(camera, event, n)

% NOTE: this function is currently overlapping with endOfTrial and saveTrial
% Need to use this version with multiple triggered cameras

gui = getappdata(0, 'gui'); 
handles = guidata(gui.maingui);

src = getselectedsource(camera);
metadata = getappdata(0,'metadata');

if get(handles.checkbox_record, 'Value') == 1  
    
end

[vid, vid_ts] = getdata(camera, camera.FramesAvailable);

videoname=sprintf('%s\\%s_cam%02d_%03d', metadata.folder, metadata.basename, n, metadata.cam(n).trialnum);

save(videoname, 'vid', 'vid_ts', '-v6')

lastvideonum = sprintf('lastvideo%d', n);
setappdata(0, lastvideonum, vid);

fprintf('Data from camera %d trial %03d successfully written to disk.\n', n, metadata.cam(n).trialnum)

if isprop(src,'FrameStartTriggerSource')
    src.FrameStartTriggerSource = 'Freerun';  % Switch from TTL to free run mode
else
    src.TriggerSource = 'Freerun';
end