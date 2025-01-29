function unfreezeCam()

cameras = getappdata(0,'cameras');

for i=1:length(cameras)
    stop(cameras(1))
end

