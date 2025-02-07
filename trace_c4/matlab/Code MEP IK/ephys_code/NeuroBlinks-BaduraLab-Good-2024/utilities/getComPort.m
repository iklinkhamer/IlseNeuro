function comport = getComPort(device_id)

% Requires pio command to be in system path
% pio is found in PlatformIO Core CLI
[st,out] = system('pio device list');

if st ~= 0
    comport = [];
end

lines = splitlines(out);

comports = regexp(lines, '^COM\d+', 'match', 'once');
devices = regexp(lines, device_id, 'match', 'once');

comport_matches = find(cellfun(@(a) ~isempty(a),comports));
device_matches = find(cellfun(@(a) ~isempty(a),devices));

idx = comport_matches(find(comport_matches < device_matches, 1, 'last'));

comport = string(comports{idx});