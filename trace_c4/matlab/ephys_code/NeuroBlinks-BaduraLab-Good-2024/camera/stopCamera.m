function stopCamera(n)

    % NEED TO GENERALIZE TO "N" CAMERAS
cameras = getappdata(0,'cameras');

% Stop video object, callback will then execute
stop(cameras{n})