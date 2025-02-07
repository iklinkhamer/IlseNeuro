function startCamera(n)

cameras = getappdata(0,'cameras');
metadata = getappdata(0,'metadata');
config = getappdata(0, 'config');

src = getselectedsource(cameras(n));

% % Need to set different property here depending on camera capabilities
% if isprop(src, 'FrameStartTriggerSource')
%     src.FrameStartTriggerSource = config.camera(n).triggermode;
% elseif isprop(src, 'TriggerSource')
%     src.TriggerSource = config.camera(n).triggermode;
% else
%     % Do nothing. Camera doesn't have this property?
%     warning('Camera %d does not have Trigger Source property', n)
% end
src.TriggerSource='Line1';
if n==1 %NHE only cam1 can be triggered
%src.TriggerMode = 'On';
end %NHE

% TODO: REMOVE THE FOLLOWING LINE BEFORE PRODUCTION
% src.TriggerMode = 'Off';    % For now override hardware trigger so we don't need Arduino connected

% cameras(n).TriggerRepeat = 0;
cameras(n).FramesPerTrigger = metadata.cam(n).frame_rate * (sum(metadata.cam(1).time) / 1e3);   % Always reference to camera 1 time

if n==1 % for Ace camera only
 src.AcquisitionBurstFrameCount = ceil(metadata.cam(n).frame_rate * (sum(metadata.cam(1).time) / 1e3));    % Note this is limited to 255 frames (see documentation for FrameBurstStart)
end

% if strcmp(cameras(n).triggermode,'Manual to disk')
%     % We have to set up the video writer object here 
%     basename = sprintf('%s\\%s_camera%02d', metadata.folder, metadata.basename, n);
%     videoname=sprintf('%s_%03d.mp4', basename, metadata.cam(n).trialnum);
%     diskLogger = VideoWriter(videoname,'MPEG-4');
%     set(diskLogger,'FrameRate',20);
%     cameras(n).DiskLogger = diskLogger;
% else
%     cameras(n).DiskLogger = [];
% end

% cameras(n).StopFcn = {@stopCameraCallback, n};
if n == 1   % Use Camera 1 as master timer of acquisition
    cameras(n).StopFcn = @endOfTrial;
end

flushdata(cameras(n));
start(cameras(n));

% Trigger Trial
%---------------- NEED TO TRIGGER TRIAL HERE --------------------%

