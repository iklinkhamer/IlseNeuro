% Rig specific settings
config.ALLOWEDDEVICES = {'arduino'};

% Better to use hash table?
% index into "camera" is camera # for particular rig, index into IDS is rig number
config.camera(1).IDS = {'02-2020C-07023'}; %{'22798789', '22804685', '22804686', '22486477'};  % DeviceInfo.DeviceName returned from imaqhwinfo
config.camera(2).IDS = {'acA1280-60gc'}; %{'22468629', '22468628', '22804762', '22784550'}; % DeviceInfo.DeviceName returned from imaqhwinfo

config.camera(1).triggermode = 'Line1';     % camera1
config.camera(2).triggermode = 'Software'; % camera2
config.camera(1).binning = 2;   % Horizontal and vertical binning for sensor (reduces resolution)
config.camera(2).binning = 2;

config.camera(1).thresh = 50/255;
config.camera(2).thresh = 50/255;

% config.camera(1).roiposition = [528 89 640 512];   % [0 0 640 512] or [0 0 1280 1024]
config.camera(1).roiposition = [0 0 328 248];   % [0 0 640 512] or [0 0 1280 1024]
% config.camera(2).roiposition = [468 8 640 480];   % [0 0 640 480] or [0 0 1280 960]
config.camera(2).roiposition = [0 0 640 480]; %!!cannot be fullsize for GigE cam, overloads bandwidth causing missing
%frames
config.camera(1).fullsize = [0 0 1280 1024];   % [0 0 640 512] or [0 0 1280 1024]
config.camera(2).fullsize = [0 0 640 480];   % [0 0 640 480] or [0 0 1280 960]

% Consider making this a per mouse setting
% config.camera(1).eyelid_roi = [401 63 67 58];
%config.camera(1).eyelid_roi = [426 149 44 48];
config.camera(1).eyelid_roi = [100 100 100 100];

config.tube_delay(1) = 12; %NHE compensates for delay in time needed for air to flow through tube
config.tube_delay(2) = 0;
config.tube_delay(3) = 0;
config.tube_delay(4) = 0;

% Index is rig number
% Arduino Zero
 config.MICROCONTROLLER_IDS = { 'USB\VID_2A03&amp;PID_003D\55431303937351E0A1D1' };% 'D2-D53', 'A0-A11','2F2E6A78', '3B3B06F1'};
% Teensy
%config.MICROCONTROLLER_IDS = {'4788080', '3580040', '4362910', '5332140'};

% For whitenoise devices, typically Teensy 3.2 with audioshield (zero disables)
config.WHITENOISE_DEVICE_IDS = {'5302470', '5404780', '5140440', '5139470'};
% config.WHITENOISE_DEVICE_IDS = {'', '', '', ''};