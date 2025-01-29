%% I.K. 6-6-24
% Called by function: ephysOverviewIK
function  MouseStruct = getBehaviorEyeDataEphys(mouseCodes, Gene, loc, Type)
P = IkUtils.getParams();
mice = defaultMice();
if nargin < 2
    if exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s.mat", mouseCodes)), 'file')
        dataStruct = load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s.mat", mouseCodes)));
        MouseStruct = dataStruct.(string(fieldnames(dataStruct)));
        return       
    end
else
    if exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s%s%s.mat", Gene, loc, Type)), 'file')
        dataStruct = load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s%s%s.mat", Gene, loc, Type)));
        MouseStruct = dataStruct.(string(fieldnames(dataStruct)));
        return
    end
end
type = ones(1,P.nTrials);
Mousetraces=nan(P.n_sessions,size(mouseCodes,2)*P.nTrials,P.nTimeSteps);
traces_per_mouse = struct;
for m = 1 : length(mouseCodes)
    mcode = mouseCodes(m);
    ephys_dates = mice([mice.code] == mcode).ephysdates;
    ephys_dates_ = split(ephys_dates,'-');
    path_list = getBehaviorPathList(mcode);
    if ~isempty(ephys_dates)
        parts = split(ephys_dates, "-"); 
                        if numel(ephys_dates) > 1
                            ephysDatesJoined = join(parts, "", 3);
                        else
                            ephysDatesJoined = join(parts, "");
                        end
                        ephysDatesJoined(ephysDatesJoined.strlength > 6) = extractAfter(ephysDatesJoined(ephysDatesJoined.strlength > 6), 2);

       
        ephys_dates_trimmed = ephysDatesJoined;
        [~, folders_] = fileparts(path_list);
        path_list_ = path_list(cellfun(@(x) ~isnan(str2double(x)), folders_));
        [~, folders] = fileparts(path_list_);
        path_list = path_list_(find(contains(folders, ephys_dates_trimmed(1))):find(contains(folders,ephys_dates_trimmed(end))));
    end

    lastSession = min(P.n_sessions, length(path_list));

    for f = 1 : length(path_list)
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mcode, s, folder);
        try
            load(fullfile(folder, filename(end)));
        catch
            continue
        end
        if length(behavior_trial_data.session_of_day) < P.nTrials
            addNantoTrialData(mcode);
            load(fullfile(folder, filename(end)));
        end
        try
            different_fields = setdiff(fieldnames(behavior_trial_data), fieldnames(trial_data));
            behavior_trial_data = rmfield(behavior_trial_data, different_fields);
        catch
        end
        try
            different_fields = setdiff(fieldnames(trial_data), fieldnames(behavior_trial_data));
            trial_data = rmfield(trial_data, different_fields);
        catch
        end
        trial_data(f) = behavior_trial_data;
        smoothedTraces = smootheTraces(trial_data(f).eyelidpos);
        traces_new(f).traces = normalizeEyelidTraces(smoothedTraces, trial_data(f));
        trial_data(f).tracesnorm = traces_new(f).traces;
    end

    for j=1:lastSession
        Mousetraces(j,(m-1)*P.nTrials+1:m*P.nTrials,:) = trial_data(j).tracesnorm(1:P.nTrials,:);
        traces_per_mouse(m).s(j).traces =  trial_data(j).tracesnorm(1:P.nTrials,:);
    end
    type(trial_data(1).c_csdur == 0) = 0;
    type(trial_data(1).c_usdur == 0) = 2;
    Mousetype = ones(length(mouseCodes), P.nTrials);
    for i = 1 : size(Mousetype,1)
        Mousetype(i, :) = type;
    end
    ts=trial_data(1).tm(1,:);    
end
MouseStruct.traces = Mousetraces;
MouseStruct.type = Mousetype;
MouseStruct.ts = ts;
MouseStruct.indvMouseTraces = traces_per_mouse;
MouseStruct.codes = mouseCodes;

if nargin < 2
    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s.mat", mcode)), "MouseStruct");
else
    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("ephysOverviewData_%s%s%s.mat", Gene, loc, Type)), "MouseStruct");
end
end
