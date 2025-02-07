function camera = addCam(serialNum, config)
% Create camera object for camera with given serial number

fprintf('Connecting to camera %s...\n', serialNum)

% Get list of configured cameras
%     foundcams = imaqhwinfo('gige');
foundcams = imaqhwinfo(config.CAMADAPTOR);

% % This is the previous version, which required opening a connection to each camera
% % This simplified newer version can be used instead as long as the camera serial number
% % is reported in DeviceInfo.DeviceName
% founddeviceids = cell2mat(foundcams.DeviceIDs); 

% if isempty(founddeviceids)
%     error('No cameras found')
% end

% %====== Getting camera ID  =====
% match=0;

% for i=1:length(founddeviceids)
%     camera = videoinput(config.CAMADAPTOR, founddeviceids(i), 'Mono8');
%     src = getselectedsource(camera);
%     if strcmp(src.DeviceSerialNumber,serialNum)
%         match=1;
%         break
%     end
%     % If we get here the current camera didn't match
%     delete(camera)
% end

% if ~match
%     error('The camera you specified (%s) could not be found',serialNum);
% end


founddevices = foundcams.DeviceInfo;

if isempty(founddevices)
    error('No cameras found')
end

match = 0 ;

for i=1:length(founddevices)
    if contains(founddevices(i).DeviceName, serialNum)
        match=1;
        camera = videoinput(config.CAMADAPTOR, founddevices(i).DeviceID, 'Mono8');
        break
    end
end
%ROIPosBefore=camera.ROIPosition
if ~match
    error('The camera you specified (%s) could not be found', serialNum);
end