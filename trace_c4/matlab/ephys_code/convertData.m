%% IK 7/9/23     load data and convert to RAWdata.bin file and also create timestamps.mat and stimstaps.mat files.
function convertData()
app = getUImain(scriptname = mfilename);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    disp("Calculating...")
    appData = evalin('base', 'appData');

    previewData = appData.previewData;
    saveConvertedDataLogical = appData.convertRawSpikingData;
    saveStimTimestampsLogical = appData.convertTimeStamps;
    mousecodes = appData.mousecodes;
    sessions = appData.sessions;

    for mcode = mousecodes
        for s = sessions
            mouseCodes = [defaultMice().code];
            mouseNames = [defaultMice().name];

            ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
            mouseNames = mouseNames(ephysMiceMask);
            mouseCodes = mouseCodes(ephysMiceMask);


            fprintf("Mouse: %s Session: %d \n", mcode, s);
            if saveConvertedDataLogical == 1
                saveConvertedData(mcode, s);
            end

            if previewData == 1
                previewData(mcode, s);
            end

            if saveStimTimestampsLogical == 1
                saveStimTimestamps(mcode, s);
            end
        end
    end
    disp("Loop end. Waiting for UI input...")
end
end

function saveConvertedData(mcode, s) % convert data
[files, folders] = getSortedFiles(mcode, s);
if isempty(files)
    disp("Channel.continuous files not found. Returning to main function.")
    return
end
if ~all(contains(files, mcode))
    reSaveWithStamp(mcode = mcode)
end

for ch = 1 : IkUtils.getParams().num_channels
    file = files(ch);
    folder = folders(ch);
    loaded_data = load_open_ephys_data(fullfile(folder, file)); %load data
    if ~exist('data','var')
        data(ch,:) = loaded_data;
    else
        data(ch,:) = loaded_data(1:length(data));
    end
end
stamp = nameDateStampFiles(mcode = mcode, s = s);
fileID = fopen(fullfile(folder, sprintf('RAWdata_%s.bin', stamp)),'w');

fwrite(fileID, data, 'double');
fclose(fileID);
end

function saveStimTimestamps(mcode, s) % make timestamps and stimstamps files
[files, folders] = getSortedFiles(mcode, s);

filepath = fullfile(folders(1), files(1));
if isempty(filepath)
    disp("channel.continuous file not found, returning to main function.")
    return
end
[~,timestamps(1,:), ~] = load_open_ephys_data(filepath); %load data

[file_events, folder_events] = recursiveFileSearch("all_channels*.events", mcode, s);
if isempty(file_events)
    disp("Events file not found, returning to main function.")
    return
end
[~,timestamps2, ~] = load_open_ephys_data(fullfile(folder_events(1), file_events(1)));

times = timestamps2 - timestamps(1,1);

stamp = nameDateStampFiles(mcode = mcode, s = s);

save(fullfile(folders(1),sprintf('stimtimes_%s.mat', stamp)),'times')

T = timestamps(:,1);
save(fullfile(folders(1), sprintf('timestamps_%s.mat', stamp)), 'T')
end

function previewData(mcode, s) % preview data
P = getParams();
[files, folders] = getSortedFiles(mcode, s);
if numel(files) ~= P.num_channels
    disp("Number of files found not equal to number of channels on the probe. Returning to main function.")
    return
end
for ch = 1 : P.num_channels
    filepath = fullfile(folders(ch), files(ch));
    if isempty(filepath)
        continue
    end
    [~, timestamps(ch,:), ~] = load_open_ephys_data(filepath); %load data
end
try
    dir_content = dir(fullfile(folders(1), 'data*.m'));
    dir_content = dir_content(~ismember({dir_content.name}, {'.', '..'}));
    data = load(dir_content.name);
catch
    disp("No saved data found. Returning to main function.")
    return
end

figure ()
title('Signal vs time(s)')
for ch = 1:P.num_channels
    subplot(8,round(P.num_channels/8),ch)
    plot(timestamps(ch,1:15000),data(ch,1:15000),'color', 'black')
    ylim([-400 400])
    title(ch)
    xlabel('time(s)')
end
end


function [sorted_files, sorted_folders] = getSortedFiles(mcode, s)
[files, folders] = recursiveFileSearch("100_*.continuous", mcode, s);
if ~isempty(files)
for f = 1:numel(files)
    [~, file, ~] = fileparts(files(f));
    parts = split(file, '_');
    channelPart(f) = cellfun(@str2double, regexp(parts(2), '\d+', 'match'));
    formatPart = regexp(parts(2), '[^\d]+', 'match');
end
else
    disp("No channel file found")
    sorted_files = [];
    sorted_folders = [];
    return
end

[~, new_idcs] = sort(channelPart);
sorted_files = files(new_idcs);
sorted_folders = folders(new_idcs);

if isempty(formatPart) && numel(sorted_files) == 32
    open_ephys_mapping = IkUtils.getParams().openEphysMapping;
    open_ephys_mapping(open_ephys_mapping > 16) = open_ephys_mapping(open_ephys_mapping > 16) - 32;
    [~, open_ephys_mapping] = sort(open_ephys_mapping);
    sorted_files = sorted_files(open_ephys_mapping);
    sorted_folders = sorted_folders(open_ephys_mapping);
end

diff_folders = unique(sorted_folders);
for f = 1 : numel(diff_folders)
    folders_type = sorted_folders(sorted_folders == diff_folders(f));
    if numel(folders_type) == IkUtils.getParams().num_channels
        folder_mask(f) = true;
    else
        folder_mask(f) = false;
    end
end
if sum(folder_mask) ~= 0 
    folder = diff_folders(find(folder_mask, 1));
    sorted_files = sorted_files(sorted_folders == folder); 
    sorted_folders = sorted_folders(sorted_folders == folder);
elseif numel(sorted_files) == 32
    sorted_files = [];
    sorted_folders = [];
    disp("Number of files found not equal to number of channels on the probe. Returning to main function.")
    return
end
end





