function endOfTrial(obj, event)
% This function is run when the camera is done collecting frames, then it calls the appropriate 
% function depending on whether or not data should be saved

gui = getappdata(0, 'gui'); 
cameras = getappdata(0, 'cameras');

handles = guidata(gui.maingui);

if get(handles.checkbox_record, 'Value') == 1  
    savetrial();
else
    nosavetrial();  
end


    src = getselectedsource(cameras(1));
    % % Set camera to freerun mode so we can preview
    % if isprop(src, 'FrameStartTriggerSource')
    %     src.FrameStartTriggerSource = 'Freerun';
    % elseif isprop(src, 'TriggerSource')
    %     src.TriggerSource = 'Freerun';
    % else
    %     % Do nothing. Camera doesn't have this property?
    %     warning('Camera %d does not have Trigger Source property', i)
    % end
    src.TriggerMode = 'Off';
    


% Do anything here that you want to start at the beginning of the ITI, such as setting new context

checkContext(handles)