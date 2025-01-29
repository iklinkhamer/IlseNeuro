function disconnectMicrocontroller(ucontroller_obj)

    % Set it up as a switch statement because we will add other options later
    switch config.device
        case 'arduino'
            disconnectArduino(ucontroller_obj);
        otherwise
            % do nothing
    end
