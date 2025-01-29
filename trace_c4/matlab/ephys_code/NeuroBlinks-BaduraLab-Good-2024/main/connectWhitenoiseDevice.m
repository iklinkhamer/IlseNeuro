function whitenoise_device = connectWhitenoiseDevice(config)

    %% -- start serial communication to Teensy ---
    disp('Connecting to white noise device...')

    % This is the faster version if you have PlatformIO command line tools installed
    % See connectTeensy for slower version
    com_port = getComPort(config.WHITENOISE_DEVICE_IDS{config.rig});

    if isempty(com_port)
        whitenoise_device = [];
        error('No white noise device found for requested rig (%d)', config.rig);  
    end

    whitenoise_device=serial(com_port, 'BaudRate', 115200);
    % whitenoise_device.InputBufferSize = 512*8;
    % whitenoise_device.DataTerminalReady='off';	% to prevent resetting teensy on connect
    whitenoise_device.DataTerminalReady='on';	    % Does this reset device for teensy? Seems to need to be set for device to respond to serial inputs
    fopen(whitenoise_device);

end

