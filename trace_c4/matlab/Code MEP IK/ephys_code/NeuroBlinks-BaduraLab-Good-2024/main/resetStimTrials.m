function resetStimTrials()
%    Note: this function is not currently used

trials=getappdata(0,'trials');

trials.stimnum=0;

setappdata(0,'trials',trials);

drawnow         % Seems necessary to update appdata before returning to calling function