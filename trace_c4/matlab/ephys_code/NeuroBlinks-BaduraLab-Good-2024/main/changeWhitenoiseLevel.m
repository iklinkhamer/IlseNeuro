function changeWhitenoiseLevel(new_level)

config = getappdata(0, 'config');
whitenoise_device = getappdata(0, 'whitenoise_device');

if ~isempty(whitenoise_device)
    fwrite(whitenoise_device, new_level, 'uint8')
end