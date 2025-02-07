function changeAmbientLEDIntensity(new_value)

    config = getappdata(0, 'config');
    microController = getappdata(0, 'microController');

    microControllerVariablesEnum;

    value = round(new_value * (2^config.PWM_RESOLUTION - 1));

    fwrite(microController, uController.AMBIENTLED, 'uint8');              % header
    fwrite(microController, value, 'int16');  % data