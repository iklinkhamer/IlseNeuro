function newFrameCallbackCam1(obj, event, himage)

gui = getappdata(0, 'gui'); 
handles = guidata(gui.maingui);
% cameras = getappdata(0, 'cameras');
metadata = getappdata(0, 'metadata');  

persistent timeOfStreamStart

if isempty(timeOfStreamStart)
    timeOfStreamStart = clock;
end

persistent timeSinceLastTrial

if isempty(timeSinceLastTrial)
    timeSinceLastTrial = clock;
end

persistent eyedata

if isempty(eyedata)
    eyedata = NaN * ones(500, 2);  
end

plt_range = -2100;

persistent eyeTrace

try

if isempty(eyeTrace)
    set(0, 'currentfigure', gui.maingui)
%     set(gui.maingui, 'CurrentAxes', handles.axes_eye)
%     cla
    eyeTrace = plot(handles.axes_eye, [plt_range 0], [1 1]*0, 'k-'); hold on
    set(handles.axes_eye, 'Color', [240 240 240]/255, 'YAxisLocation', 'right');
    set(handles.axes_eye, 'XLim', [plt_range 0], 'YLim', [-0.1 1.1])
    set(handles.axes_eye, 'XTick', [-3000:500:0], 'box', 'off')
    set(handles.axes_eye, 'YTick', [0:0.5:1], 'YTicklabel', {'0' '' '1'})
    set(handles.axes_eye, 'FontSize', 8);
end


% --- eye trace ---
wholeframe = event.Data;
roi = wholeframe .* uint8(metadata.cam(1).mask);
eyelidpos = sum(roi(:)>= 255 * metadata.cam(1).thresh);

% --- eye trace buffer ---
eyedata(1:end-1, :) = eyedata(2:end, :);
timeSinceStreamStartMS = round(1000 * etime(clock, timeOfStreamStart));
eyedata(end, 1) = timeSinceStreamStartMS;
eyedata(end, 2) = (eyelidpos - metadata.cam(1).calib_offset) / metadata.cam(1).calib_scale; % eyelid pos

set(eyeTrace, 'XData', eyedata(:, 1) - timeSinceStreamStartMS, 'YData', eyedata(:, 2))
set(himage, 'CData', event.Data)
        
        
% --- Check if new trial should be triggered ----

if get(handles.toggle_continuous, 'Value')  == 1
        
    stopTrial = str2double(get(handles.edit_StopAfterTrial, 'String'));
    if stopTrial > 0 && metadata.cam(1).trialnum > stopTrial
        set(handles.toggle_continuous, 'Value', 0);
        set(handles.toggle_continuous, 'String', 'Start Continuous');
    end
        
    elapsedTimeSinceLastTrial = etime(clock, timeSinceLastTrial);
    timeLeft = metadata.stim.c.ITI - elapsedTimeSinceLastTrial;
    interTrialTime=metadata.stim.c.ITI;
    set(handles.trialtimecounter, 'String', num2str(round(timeLeft)))
    
    if timeLeft <= 0
        eyeok = checkeye(handles, eyedata);
        if eyeok
            startTrial(handles)
            timeSinceLastTrial = clock;
        end
    end
end

catch ME
    disp('error in starting the trial')
    throw(ME)

end

end

function eyeok = checkeye(handles, eyedata)

    eyethrok = (eyedata(end, 2) < str2double(get(handles.edit_eyelidthresh, 'String')));
    eyedata(:, 1) = eyedata(:, 1) - eyedata(end, 1);  
    recenteye = eyedata(eyedata(:, 1) > -1000 * str2double(get(handles.edit_stabletime, 'String')), 2);
    eyestableok = ((max(recenteye) - min(recenteye)) < str2double(get(handles.edit_stableeye, 'String')));
    eyeok =  eyestableok; %eyethrok &&
   %eyeok = 1; %NHE added to circumvent stability issues due to light fluctuations during testing, delete after!
end
