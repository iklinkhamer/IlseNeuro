function ResetCamTrials()

metadata=getappdata(0,'metadata');
cameras=getappdata(0,'cameras');

for i=1:length(cameras)
    metadata.cam(i).trialnum=1;
end

setappdata(0,'metadata',metadata);

drawnow         % Seems necessary to update appdata before returning to calling function