% Configuration settings that might be different for different users/computers
% Should be somewhere in path
DEFAULTDEVICE= 'arduino';
DEFAULTRIG=1;
config.PWM_RESOLUTION = 12;        % Bits of resolution of PWM signal in microcontroller (for LED intensity control)

config.CAMADAPTOR = 'gentl';

% --- camera settings ----
% Need to rework these configuration settings for consistency in naming/numbering different cameras
config.camera(1).initExposureTime=300;  % Exposure times for the two cameras (e.g., eyelid, pupil)
config.camera(2).initExposureTime=450;  % Exposure times for the two cameras (e.g., eyelid, pupil)
config.camera(1).FrameRate = 200;   % Frame rates for the two cameras (e.g., eyelid, pupil)
config.camera(2).FrameRate = 60;   % Frame rates for the two cameras (e.g., eyelid, pupil)

% Supress know warnings
% This one is about Camera2, which only has a software trigger on the pulse
warning('off','imaq:gentl:immedOrManTriggerTurningTriggerModeOff')
 
