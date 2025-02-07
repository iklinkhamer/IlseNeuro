configureGige(camera, config)

src = getselectedsource(camera);

src.ExposureTimeAbs = config.cam.initExposureTime;
src.AllGainRaw=12;				% Tweak this based on IR light illumination (lower values preferred due to less noise)
% src.StreamBytesPerSecond=124e6; % Set based on AVT's suggestion
src.StreamBytesPerSecond=80e6; % Set based on AVT's suggestion

% src.PacketSize = 9014;		% Use Jumbo packets (ethernet card must support them) -- apparently not supported in VIMBA
src.PacketSize = 8228;		% Use Jumbo packets (ethernet card must support them) -- apparently not supported in VIMBA
src.PacketDelay = 2000;		% Calculated based on frame rate and image size using Mathworks helper function
camera.LoggingMode = 'memory'; 
src.AcquisitionFrameRateAbs=200;
camera.FramesPerTrigger=ceil(recdur/(1000/200));

% triggerconfig(camera, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
% set(src,'AcquisitionStartTriggerMode','on')
% set(src,'FrameStartTriggerSource','Freerun')
% set(src,'AcquisitionStartTriggerActivation','RisingEdge')
% set(src,'AcquisitionStartTriggerSource','Line1')

triggerconfig(camera, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
src.FrameStartTriggerMode = 'On';
src.FrameStartTriggerActivation = 'LevelHigh';

% This needs to be toggled to switch between preview and acquisition mode
% It is changed to 'Line1' in MainWindow just before triggering Arduino and then
% back to 'Freerun' in 'endOfTrial' function
% src.FrameStartTriggerSource = 'Line1';
src.FrameStartTriggerSource = 'Freerun';