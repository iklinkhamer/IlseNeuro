%% IK 26-4-24 Retrieve ephys data struct
function data = getData(mcode)
mice = defaultMice();
dates = mice([mice.code] == mcode).ephysdates;
if ~isempty(dates)
    numdates = numel(dates);
    sessions = 1:numdates;
elseif contains(mcode, "Shank2KO")
    sessions = 1:IkUtils.getParams().n_sessions;
else
    disp("No ephys dates found in defaultMice")
    data = [];
    return
end

for s = sessions
    [filename, folder] = recursiveFileSearch("StructEphysData*.mat", mcode, s);
    if numel(filename) > 1
        stamp = nameDateStampFiles(mcode = mcode, s = s);
        for f = 1 : numel(filename)
            if contains(filename(f), stamp)
                filename = filename(f);
                folder = folder(f);
                break
            end
        end
    end

    if isempty(filename)
        %fprintf("Ephys data struct not found for session: %d\n", s)
        continue
    end

    try
        load(fullfile(folder, filename))
        if ~isempty(unit) 
            if ~ isempty(unit.neuron)
                data(s) = unit;
            end
        end
    catch
        fprintf("Data file not found for session %d", s)
        if ~exist("data", 'var')
            fields_data = fields(data);
            if ~isempty(fields_data)
                for f = 1 : numel(fields_data)
                    data(s).(fields_data(f)) = [];
                end
            end
        end
    end
end
if ~exist("data", 'var')
    data = [];
end
end