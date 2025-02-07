function savetrial()

% For encoder
% 4096 counts = 1 revolution = 15.24*pi cm for 6 inch diameter cylinder
counts2cm = @(count) double(count) ./ 4096 .* 15.24 .* pi; 

microControllerVariablesEnum;

% Load objects from root app data
cameras=getappdata(0, 'cameras');
metadata = getappdata(0, 'metadata');
config = getappdata(0, 'config');

setappdata(0, 'lastmetadata', metadata);

%%
%% ------------- TODO: Save encoder data separately from camera data ---------------------- %%
%%
% Get encoder data from microController
if config.use_encoder == 1

  if isappdata(0, 'microController')
    microController = getappdata(0, 'microController');

    if microController.BytesAvailable > 0
      fread(microController, microController.BytesAvailable); % Clear input buffer
    end

    fwrite(microController, uController.REQUESTDATA, 'uint8');  % Tell microController we're ready for it to send the data

    data_header = (fread(microController, 1, 'uint8'));
    if data_header == uController.ENCODERCOUNTS
      encoder.counts = (fread(microController, 200, 'int32'));
      encoder.displacement = counts2cm(encoder.counts - encoder.counts(1));
    end

    time_header = (fread(microController, 1, 'uint8'));
    if time_header == uController.ENCODERTIME
      encoder.time = (fread(microController, 200, 'uint32'));
      encoder.time = encoder.time - encoder.time(1);
    end

    encoderdatafname = sprintf('%s\\%s_%03d_encoder', metadata.folder, metadata.basename, metadata.cam(1).trialnum);
    save(encoderdatafname, 'encoder', '-v6')

    if config.verbose
      fprintf('Encoder data for trial %03d successfully written to disk.\n', metadata.cam(1).trialnum)
    end

  end

end

% Get videos from cameras
for i=1:1 %length(cameras)

  if strcmp(cameras(i).Running, 'on')
    stop(cameras(i))
  end

  acquiredFrameDiff = cameras(i).FramesAvailable < cameras(i).FramesPerTrigger * (cameras(i).TriggerRepeat + 1);
  if acquiredFrameDiff > 0
      warning('There is a difference of %d between the number of frames you expected and the number that were acquired', acquiredFrameDiff)
  end

  [vid, vid_ts] = getdata(cameras(i), cameras(i).FramesAvailable);
  vid_ts = vid_ts - vid_ts(1);

  % Keep data from last trial in memory even if we don't save it to disk
  setappdata(0, sprintf('lastvid%d', i), vid);

  videoname = sprintf('%s\\%s_%03d_cam%d', metadata.folder, metadata.basename, metadata.cam(i).trialnum, i);
  % if exist('encoder', 'var')
  %     save(videoname, 'vid', 'vid_ts', 'metadata', 'encoder', '-v6')
  % else
  %     save(videoname, 'vid', 'vid_ts', 'metadata', '-v6')
  % end
  save(videoname, 'vid', 'vid_ts', 'metadata', '-v6')

  if config.verbose
    fprintf('Video data from camera %d for trial %03d successfully written to disk.\n', i, metadata.cam(i).trialnum)
  end

  if i == 1
    onlineBehaviorAnalysis(vid, metadata);
  end

  metadata.cam(i).trialnum = metadata.cam(i).trialnum + 1;

end

% online_bhvana(vid);

setappdata(0, 'metadata', metadata);

drawnow         % Seems necessary to update appdata before returning to calling function
