function camera = configureBaslerPulse(camera, config)

src = getselectedsource(camera);

src.AcquisitionFrameRate = config.camera(2).FrameRate;   % index 1 is for high speed (eyelid) cam [this is hacky and should change in later version]

src.ExposureAuto = 'off';
src.ExposureTime = config.camera(2).initExposureTime;

src.GainAuto = 'off';
src.Gain=0;				% Tweak this based on IR light illumination (lower values preferred due to less noise)
src.Gamma = 1.5;

src.BinningHorizontal=config.camera(2).binning;
src.BinningVertical=config.camera(2).binning;

camera.ROIPosition = config.camera(2).roiposition;

camera.LoggingMode = 'memory'; 
camera.FramesPerTrigger = ceil(config.trial_length_ms / (1000 / config.camera(2).FrameRate));

triggerconfig(camera, 'Immediate', 'none', 'none');
src.TriggerMode = 'Off';
src.TriggerActivation = 'RisingEdge';

src.TriggerSource = 'Software';