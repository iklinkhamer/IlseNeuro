function camera = initCam(config)

% First delete any existing image acquisition objects
imaqreset

disp('Finding camera...')
 % Get list of configured cameras
foundcams = imaqhwinfo(config.CAMADAPTOR);
founddeviceids = cell2mat(foundcams.DeviceIDs); 

if isempty(founddeviceids)
    error('No cameras found')
end

%====== Getting camera ID  =====
cam_num = 0;
%     cam_num = 1;

for i=1:length(founddeviceids)
    camera = videoinput(config.CAMADAPTOR, founddeviceids(i), 'Mono8');
    src = getselectedsource(camera);
    if strcmp(src.DeviceID, config.camera(1).IDS{config.rig})
        cam_num = i;
    end
    delete(camera)
end

if ~cam_num
    error('The camera you specified (%d) could not be found', config.rig);
end

disp('Creating camera object...')

camera = videoinput(config.CAMADAPTOR, cam_num, 'Mono8');

configureBaslerAce(camera, config)

% setappdata(0,'camera',camera)