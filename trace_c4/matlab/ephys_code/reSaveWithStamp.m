%% IK 30-04-2024
function reSaveWithStamp(kwargs)
arguments
    kwargs.directory = '/home/i.klinkhamer/Documents/Data/spikeSortingUnits/';
    kwargs.mcode = '';
    kwargs.allSessions = true
end
directory = kwargs.directory;
P = IkUtils.getParams();

if isempty(kwargs.mcode)
    app = getUImain(scriptname=mfilename);
end

while true
    if isempty(kwargs.mcode)
        uiwait(app.UIFigure)
        if ~isvalid(app)
            % If the UI figure has been closed, exit the loop
            break;
        end
        appData = evalin('base', 'appData');
        mousecodes = appData.mousecodes;
        sessions = appData.sessions;
    else
        mousecodes = kwargs.mcode;
        if kwargs.allSessions == 0
            prompt2 = "Enter session number: ";
            sessions = input(prompt2);
        else
            mice = defaultMice();
            dates = mice([mice.code] == kwargs.mcode).ephysdates;
            numdates = sum(~cellfun(@isempty, dates));
            sessions = 1:numdates;
        end

    end

    % Loop through each file and rename it
    for mcode = mousecodes
        files = dir(fullfile(directory, mcode));
        if ~isempty(files)
            renameMoveFiles(files, mcode)
        end
        for s = sessions
            files = dir(fullfile(directory, mcode, P.s(s)));
            renameMoveFiles(files, mcode, s)
        end
    end
    
    if ~isempty(kwargs.mcode)
        break;
    else
        disp("End of loop. Waiting for new UI input.")
    end
end
end

function renameMoveFiles(files, mcode, s)
if isempty(files)
    return
end
if nargin < 3
    stamp = nameDateStampFiles(mcode = mcode);
else
    stamp = nameDateStampFiles(mcode = mcode, s = s);
end

for i = 1:numel(files)
    if ~files(i).isdir
        % Split filename and extension
        [~, fileName, fileExt] = fileparts(files(i).name);
        fileFolder = files(i).folder;
        if contains(fileName, stamp)
            continue
        end
        parts = split(fileName, '_');
        if any(strcmp(parts, "Config"))
            idx = 1;
            for p = 1:numel(parts)
                if ~isempty(parts{p})
                    partsNoEmptyEntries{idx} = parts{p};
                    idx = idx+1;
                end
            end
            place = find(~cellfun(@isempty,regexp(partsNoEmptyEntries, "dat")) | ~cellfun(@isempty,regexp(partsNoEmptyEntries, "bin")));
            if isempty(place)
                place = find(~cellfun(@isempty,regexp(partsNoEmptyEntries, "oe")));
            end
            file_pattern = strjoin(partsNoEmptyEntries(1:place), '_');
            file_pattern_end = strjoin(partsNoEmptyEntries(place+1:end), '_');
            if ~isempty(file_pattern_end)
                file_end = "_" + string(file_pattern_end) + fileExt;
            else
                file_end = fileExt;
            end
        end
        if any(strcmp(parts, mcode)) && ~any(strcmp(parts, "Config"))
            continue
        end

        % Construct the new filename with the new part inserted between original filename and extension
        if ~any(strcmp(parts, "Config"))
            file_pattern = fileName;
            file_end = fileExt;
        end
        if nargin < 3
            newFileNameStamp = nameDateStampFiles(mcode = mcode, file_pattern = file_pattern, file_extension = file_end);
        else
            newFileNameStamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = file_pattern, file_extension = file_end);
        end
        newFileName = fullfile(fileFolder, newFileNameStamp);
        % Move the file to the new filename
        movefile(fullfile(fileFolder, files(i).name), newFileName);
    end
end

end
