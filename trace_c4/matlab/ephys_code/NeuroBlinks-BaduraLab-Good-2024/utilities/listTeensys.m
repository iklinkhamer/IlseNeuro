function arduino_list = listArduinos()
    
    % Call external function 'wmicGet' to pull in PnP info
    infostruct = wmicGet('Win32_PnPEntity');
    
    names={infostruct.Name};  % roll struct field into cell array for easy searching
    device_ids={infostruct.DeviceID};  % roll struct field into cell array for easy searching
    
    idx = find(contains(names,'USB Serial Device'));   % All devices with "USB Serial Device" in the name field (contains was introduced in Matlab recently)
    
    if isempty(idx)
        return
    end
    
    arduino_names = names(idx);
    arduino_device_ids = device_ids(idx);

    arduino_list = {};
    for i=1:length(idx)
        arduino_list{i} = sprintf('%s: %s', arduino_names{i}, arduino_device_ids{i});
    end

end