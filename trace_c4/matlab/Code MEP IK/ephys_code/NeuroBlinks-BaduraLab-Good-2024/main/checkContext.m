function checkContext(handles)
   
metadata = getappdata(0, 'metadata');
config = getappdata(0, 'config');

% If we're doing per trial white noise control
if ~strcmp(config.WHITENOISE_DEVICE_IDS{config.rig}, '') && strcmp(handles.popupmenu_whitenoise_control.String{handles.popupmenu_whitenoise_control.Value}, 'Trial table')
    
    trialvars = readTrialTable(metadata.cam(1).trialnum);
    metadata.stim.c.whitenoise = trialvars(7);

    % Tell whitenoise device new level - device only toggles if level is different
    changeWhitenoiseLevel(metadata.stim.c.whitenoise);

else
    
    metadata.stim.c.whitenoise = -1;   % -1 to indicate that whitenoise isn't being controlled by trial table (or should we use NaN?)
    
end

setappdata(0, 'metadata', metadata);