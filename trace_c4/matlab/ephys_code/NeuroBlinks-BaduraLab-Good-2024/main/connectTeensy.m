function [teensy, config] = connectTeensy(config)

    %% -- start serial communication to arduino ---
    disp('Connecting to Teensy board...')

    % % This is the slow version if you don't have PlatformIO command line tools installed (uses included wmicGet command)
    % com_ports = findTeensys(config.MICROCONTROLLER_IDS);

    % if isempty(com_ports{config.rig})
    %     error('No Arduino found for requested rig (%d)', config.rig);
    % end

    % config.com_ports(config.rig) = com_ports{config.rig};

    % This is the faster version if you have PlatformIO command line tools installed
    com_port = getComPort(config.MICROCONTROLLER_IDS{config.rig});

    if ~isempty(com_port)
        config.com_ports(config.rig) = com_port; 
    else
        error('No Teensy found for requested rig (%d)', config.rig);        
    end

    teensy=serial(config.com_ports(config.rig), 'BaudRate', 115200);
    teensy.InputBufferSize = 512*8;
    % teensy.DataTerminalReady='off';	% to prevent resetting teensy on connect
    teensy.DataTerminalReady='on';	    % Does this reset device for teensy Zero? Seems to need to be set for device to respond to serial inputs
    fopen(teensy);

end

function com_ports = findTeensys(ids)
    com_ports = cell(size(ids));
    
    % Call external function 'wmicGet' to pull in PnP info
    infostruct = wmicGet('Win32_PnPEntity');
    
    names={infostruct.Name};  % roll struct field into cell array for easy searching
    device_ids={infostruct.DeviceID};  % roll struct field into cell array for easy searching
    
    idx = find(contains(names, 'USB Serial Device'));   % All devices with "Atmel" in the name field (contains was introduced in Matlab recently)
    
    if isempty(idx)
        return
    end
    
    teensy_names = names(idx);
    teensy_device_ids = device_ids(idx);
    
    for i=1:length(ids)
        % Figure out which line the ID appears on...
        match = strfind(teensy_device_ids, ids{i});
        idx = find(~cellfun(@isempty, match));
        if ~isempty(idx)
            %...and find the corresponding COM port on that line
            com_ports{i} = regexp(teensy_names(idx), '(COM\d+)', 'match', 'once');
        end
    end
    
end