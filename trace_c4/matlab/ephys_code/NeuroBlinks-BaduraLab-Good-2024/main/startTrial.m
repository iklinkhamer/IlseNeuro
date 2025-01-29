function startTrial(handles)

metadata = getappdata(0, 'metadata');
timestamp = getappdata(0, 'timestamp');
cameras = getappdata(0, 'cameras');
microControllerVariablesEnum;

% pause(0)
% cameras(2).DiskLoggerFrameCount
% cameras(2).FramesAcquired



metadata = refreshParams(handles, metadata);
ok = uploadParams(metadata);

if ~ok
    msgbox('Problem communicating with Microcontroller, aborting trial')
    return
end

metadata.ts(2) = etime(clock, datevec(metadata.ts(1))) % IK uncomment. 
%disp('check1')
%pause(0.1)
timestamp.trial_time(timestamp.counter) = cameras(2).FramesAcquired; %etime(clock, datevec(metadata.ts(1)));
timestamp.trial_time_seconds(timestamp.counter) = etime(clock, datevec(metadata.ts(1)));
timestamp.counter = timestamp.counter +1;

startCamera(1);
%startCamera(2); %turned off as we want continuous recording
metadata.ts(2) = etime(clock, datevec(metadata.ts(1))) % IK added 1/3/24
timestamp.trial_time_seconds_after_init_cam_1(timestamp.counter) = metadata.ts(2);% IK added 1/3/24
timestamp.trial_time_frames_after_init_cam_1(timestamp.counter) = cameras(2).FramesAcquired;% IK added 1/3/24

% --- trigger via microController --
microController = getappdata(0, 'microController');
fwrite(microController, uController.RUN, 'uint8');

timestamp.trial_time_seconds_after_trigger(timestamp.counter) = metadata.ts(2); % IK added 1/3/24
timestamp.trial_time_frames_after_trigger(timestamp.counter) = cameras(2).FramesAcquired; %IK added 1/3/24
% ---- write status bar ----
% TODO: will need to rewrite this when trial vars start using enum
set(handles.text_status, 'String', sprintf('Total trials: %d\n', metadata.cam(1).trialnum));
if strcmpi(metadata.stim.type, 'conditioning')
    trialvars = readTrialTable(metadata.cam(1).trialnum + 1);
    csdur = trialvars(1);
    csnum = trialvars(2);
    isi = trialvars(4);
    usdur = trialvars(5);
    usnum = trialvars(6);
    cstone = str2num(get(handles.edit_tone, 'String'));
    if length(cstone) < 2, cstone(2) = 0; end
    
    str2 = [];
    if ismember(csnum, [5 6])
        str2 = [' (' num2str(cstone(csnum - 4)) ' KHz)'];
    end
        
    str1 = sprintf('Next:  No %d,  CS ch %d%s,  ISI %d,  US %d, US ch %d', metadata.cam(1).trialnum + 1, csnum, str2, isi, usdur, usnum);
    set(handles.text_disp_cond, 'String', str1)
end

setappdata(0, 'metadata', metadata);
setappdata(0, 'timestamp' , timestamp);
drawnow         % Seems necessary to update appdata before returning to calling function
