function onlineBehaviorAnalysis(vid, metadata)
% Eyelid trace is saved to memory (trials and metadata) for online analysis

trials=getappdata(0,'trials');
if isfield(trials,'eye') 
    if length(trials.eye) > metadata.cam(1).trialnum + 1
        trials.eye=[]; 
    end
end

if metadata.cam(1).trialnum == 1
    trials.eye=[]; 
end

% ------ eyelid trace, which will be saved to 'trials' ---- 
[trace,time]=vid2eyetrace(vid,metadata,metadata.cam(1).thresh);
trace=(trace-metadata.cam(1).calib_offset)/metadata.cam(1).calib_scale;

trials.eye(metadata.cam(1).trialnum).time=time*1e3;
trials.eye(metadata.cam(1).trialnum).trace=trace;
trials.eye(metadata.cam(1).trialnum).stimtype=lower(metadata.stim.type);
trials.eye(metadata.cam(1).trialnum).isi=NaN;


trials.eye(metadata.cam(1).trialnum).stimtime.st{1}=0; % for CS
trials.eye(metadata.cam(1).trialnum).stimtime.en{1}=metadata.stim.c.csdur;
trials.eye(metadata.cam(1).trialnum).stimtime.cchan(1)=3;
trials.eye(metadata.cam(1).trialnum).stimtime.st{2}=metadata.stim.c.isi;
trials.eye(metadata.cam(1).trialnum).stimtime.en{2}=metadata.stim.c.usdur; % for US
trials.eye(metadata.cam(1).trialnum).stimtime.cchan(2)=2;
trials.eye(metadata.cam(1).trialnum).isi=metadata.stim.c.isi;
    

% % --- this may be useful for offline analysis ----
% metadata.eye.ts0=time(1)*1e3;
% metadata.eye.ts_interval=mode(diff(time*1e3));
% metadata.eye.trace=trace;

% --- save results to memory ----
setappdata(0,'trials',trials);

drawnow         % Seems necessary to update appdata before returning to calling function


