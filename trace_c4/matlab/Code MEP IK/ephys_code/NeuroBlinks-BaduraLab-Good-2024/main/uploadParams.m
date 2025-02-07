function ok = uploadParams(metadata)

config = getappdata(0, 'config');

ok = 0; % Flag changed to 1 when function successfully completes

% Load enumerated variable indices into scope
% Enum is stored in struct
microControllerVariablesEnum;

% TODO: Need to refactor this part to make it easier to maintain
% Use same enumeration scheme as "refreshParams"
dataBuffer = zeros(1, length(fieldnames(uController)), 'int16');   % Get fieldnames for number of enums because stored in struct

dataBuffer(uController.PRETIME) = metadata.cam(1).time(1);
dataBuffer(uController.POSTTIME) = metadata.cam(1).time(2);
if strcmpi(metadata.stim.type, 'puff')
    dataBuffer(uController.USDUR) = metadata.stim.p.puffdur;
    dataBuffer(uController.USNUM) = 3;    % This is the puff channel
elseif  strcmpi(metadata.stim.type, 'conditioning')
    dataBuffer(uController.CSNUM) = metadata.stim.c.csnum;
    dataBuffer(uController.CSDUR) = metadata.stim.c.csdur;
    dataBuffer(uController.CSINT) = round(metadata.stim.c.csintensity * (2^config.PWM_RESOLUTION - 1));     % Scale to full analogWrite units
    dataBuffer(uController.USDUR) = metadata.stim.c.usdur;
    dataBuffer(uController.ISI) = metadata.stim.c.isi - config.tube_delay(config.rig);   
    if ismember(metadata.stim.c.csnum, [5, 6])
        dataBuffer(uController.TONEFREQ) = metadata.stim.c.cstone(metadata.stim.c.csnum - 4);
    end
    if ismember(metadata.stim.c.usnum, [5, 6])
        dataBuffer(uController.TONEFREQ) = metadata.stim.c.cstone(metadata.stim.c.usnum - 4);
    end
    dataBuffer(uController.USNUM) = metadata.stim.c.usnum;
    % dataBuffer(11) = metadata.stim.l.delay;
    % dataBuffer(12) = metadata.stim.l.dur;
    % dataBuffer(13) = metadata.stim.l.amp;
end

% ---- send data to microController ----
% TODO: Refactor micro controller communication to use handshake, etc
microController = getappdata(0, 'microController');
for i = 1:length(dataBuffer)
   % try
        % Note that the following code as written is inefficient because it writes a lot of empty data--need to only cycle through dataBuffer indices with values
        fwrite(microController, i, 'uint8');              % header
        fwrite(microController, dataBuffer(i), 'int16');  % data
%     catch
%         warning('Problem writing data to microcontroller: Index=%d, Value=%d\n', i, dataBuffer(i))
%         return
%     end
    % if mod(i, 4) == 0
    %     pause(0.010);
    % end
end

ok = 1;