function camera = configureBaslerAce2(camera, config)

src = getselectedsource(camera);

src.AcquisitionFrameRateAbs = config.camera(2).FrameRate;  % index 1 is for high speed (eyelid) cam [this is hacky and should change in later version]
src.AcquisitionFrameRateEnable = 'True'; %'True';
src.ExposureAuto = 'off';
src.ExposureTimeAbs = config.camera(2).initExposureTime;
src.GainAuto = 'off';
src.GainRaw=0;	
src.TriggerSource = 'Software';


camera.ROIPosition = config.camera(2).roiposition;

camera.LoggingMode = 'disk'; 
camera.TriggerRepeat = Inf; %NHE
camera.FramesPerTrigger = Inf;
%src.AcquisitionFrameRateAbs = 30;

%src.ExposureAuto = 'Continuous';
%src.ExposureTimeAbs = config.camera(2).initExposureTime;

%src.SensorReadoutMode = 'Fast';

			% Tweak this based on IR light illumination (lower values preferred due to less noise)
%src.Gamma = 1.5;

% src.BinningHorizontal=config.camera(1).binning;
% src.BinningVertical=config.camera(1).binning;
% src.BinningHorizontalMode='Sum';    % 'Sum' or 'Average'
% src.BinningVerticalMode='Sum';


%camera.FramesPerTrigger = 120; %60; %ceil(config.trial_length_ms / (1000 / config.camera(1).FrameRate))

% triggerconfig(camera, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
% set(src,'AcquisitionStartTriggerMode','on')
% set(src,'FrameStartTriggerSource','Freerun')
% set(src,'AcquisitionStartTriggerActivation','RisingEdge')
% set(src,'AcquisitionStartTriggerSource','Line1')

% triggerconfig(camera, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
% src.TriggerMode = 'Off';
% src.TriggerActivation = 'RisingEdge';
% src.TriggerSelector = 'FrameStart';
%propinfo(src, 'AcquisitionFrameCount');
%src.AcquisitionFrameCount = 50; %ceil(config.trial_length_ms / (1000 / config.camera(1).FrameRate));    % Note this is limited to 255 frames (see documentation for FrameBurstStart)


% This needs to be toggled to switch between preview and acquisition mode
% It is changed to 'Line1' in MainWindow just before triggering Arduino and then
% back to 'Freerun' in 'endOfTrial' function
% src.FrameStartTriggerSource = 'Line1';

