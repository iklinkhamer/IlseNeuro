ii = 1;
for i = 30
deviceID = imaqhwinfo('gentl');
camera = videoinput('gentl', 1,  'Mono8');
camera.LoggingMode = 'disk'; 
camera.TriggerRepeat = Inf;
camera.FramesPerTrigger = Inf;

ctime = clock;
videoname=sprintf('%s\\%s_%04d.mp4','I:\Nynke\Neuroblinks-Baduralab', 'test_cam',round(ctime(6)))
logfile = VideoWriter(videoname, 'MPEG-4');

camera.DiskLogger = logfile;

src = getselectedsource(camera);

src.AcquisitionFrameRateAbs = i;
src.AcquisitionFrameRateEnable = 'True';


src.ExposureAuto = 'off';
src.ExposureTimeAbs = 500;

%src.TriggerSelector = 'AcquisitionStart';

src.TriggerSource = 'Software';


ResultingFrameRate = src.ResultingFrameRateAbs;
start(camera);
pause(60);
delete(camera);

disp('done')

v = VideoReader(videoname);
video = read(v);
numFrames(ii) = size(video,4);
ii = ii+1;
end