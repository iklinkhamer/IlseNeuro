function nosavetrial()
% We set a callback function and trigger camera anyway so we can get instant replay

cameras = getappdata(0, 'cameras');
% pause(1e-3)
metadata = getappdata(0, 'metadata');
setappdata(0, 'lastmetadata', metadata);

for i=1:length(cameras)

    if strcmp(cameras(i).Running, 'on')
        stop(cameras(i))
    end

    acquiredFrameDiff = cameras(i).FramesAvailable < cameras(i).FramesPerTrigger * (cameras(i).TriggerRepeat + 1);
    if acquiredFrameDiff > 0
        warning('There is a difference of %d between the number of frames you expected and the number that were acquired', acquiredFrameDiff)
    end

    vid = getdata(cameras(i), cameras(i).FramesAvailable);

    % Keep data from last trial in memory even if we don't save it to disk
    setappdata(0, sprintf('lastvid%d', i), vid);

end

fprintf('Data from last trial saved to memory for review.\n')