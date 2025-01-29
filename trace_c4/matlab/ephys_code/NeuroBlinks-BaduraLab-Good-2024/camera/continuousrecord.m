function continuousrecord(handles)
n=2;
cameras = getappdata(0,'cameras');
src = getappdata(0,'src');
metadata = getappdata(0,'metadata');
%timestamp = getappdata(0, 'timestamp'); % IK added 1/3/24

% % Switch from memory logging to disk logging
% set(cameras(n), 'LoggingMode', 'disk')

%check settings
metadata.folder;
metadata.basename;
% Create a VideoWriter object with the profile set to MPEG-4
ctime = etime(clock,datevec(metadata.ts(1)));
videoname=sprintf('%s\\%s_%05d.mp4',metadata.folder,metadata.basename,round(ctime));
logfile = VideoWriter(videoname, 'MPEG-4');
timestamp.timeRightBeforeInitiationCam2 = ctime;
%tic
% Configure the video input object to use the VideoWriter object.
cameras(2).DiskLogger = logfile;

% Now that the video input object is configured forlogging data to a Motion JPEG 2000 file, initiate the acquisition.
start(cameras(2))

%toc
timestamp.acquisitionInitiationTimeCam2 = etime(clock,datevec(metadata.ts(1))); % IK added 1/3/24
% name = sprintf('%s\\%s_%s.mat',metadata.folder, metadata.basename, 'acquisitionInitiationTimeCam2')
%    save(name, 'atime');


% pause(1)
% cameras(2).DiskLoggerFrameCount
% cameras(2).FramesAcquired
timestamp.trial_time(1)=0;
timestamp.counter=1;
setappdata(0, 'timestamp', timestamp)

% % Wait for the acquisition to finish.
% wait(cameras(2), 5)
% 
% % When logging large amounts of data to disk, disk writingoccasionally lags behind the acquisition. To determine whether allframes are written to disk, you can optionally use the DiskLoggerFrameCount property.
% while (cameras(2).FramesAcquired ~= cameras(2).DiskLoggerFrameCount) 
%     pause(.1)
% end
% 
% % You can verify that the FramesAcquired and DiskLoggerFrameCount propertieshave identical values by using these commands and comparing the output.
% cameras(2).FramesAcquired
% cameras(2).DiskLoggerFrameCount 
% 
% % When the video input object is no longer needed, deleteit and clear it from the workspace.
% delete(cameras)
% disp('cam 2 deleted')
% clear camera