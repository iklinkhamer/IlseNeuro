function [ucontroller_obj, config] = connectMicrocontroller(config)

    % Set it up as a switch statement because we will add other options later
    switch config.device
        case 'arduino'
            [ucontroller_obj, config] = connectArduino(config);
        case 'teensy'
            [ucontroller_obj, config] = connectTeensy(config);
        otherwise
            ucontroller_obj = [];
    end

