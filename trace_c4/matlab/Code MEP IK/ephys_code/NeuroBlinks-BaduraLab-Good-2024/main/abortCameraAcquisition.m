function abortCameraAcquisition()

cameras = getappdata(0, 'cameras');

for i=1:length(cameras)
    src = getselectedsource(cameras(i));

    stop(cameras(i));
    flushdata(cameras(i));

    % src.FrameStartTriggerSource = 'Freerun';
    src.TriggerMode = 'off';
end