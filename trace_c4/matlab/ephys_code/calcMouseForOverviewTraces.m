%% I.K. 1-6-24
function  MouseStruct = calcMouseForOverviewTraces(mouseCodes, Gene, loc, Type)
P = IkUtils.getParams();
mice = defaultMice();
if nargin < 2
    if exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s.mat", mouseCodes)), 'file')
        dataStruct = load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s.mat", mouseCodes)));
        MouseStruct = dataStruct.(string(fieldnames(dataStruct)));
        if isfield(MouseStruct, "crampSEM")
            return
        end
    end
else
    if exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)), 'file')
        dataStruct = load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)));
        MouseStruct = dataStruct.(string(fieldnames(dataStruct)));
        if isfield(MouseStruct, "crampSEM")
            return
        end
    end
end

type = ones(1,P.nTrials);
Mousetraces=nan(P.lastDay,size(mouseCodes,2)*P.nTrials,P.nTimeSteps);
MouseCRamp=nan(P.lastDay,size(mouseCodes,2),P.nCSplusRegTrials);
Mousecramp5= nan(P.lastDay,size(mouseCodes,2), P.nCSplusRegTrials);
Mousemean=nan(P.lastDay,P.nTimeSteps);
Mouseonsets=nan(P.lastDay,size(mouseCodes,2)*P.nCSplusRegTrials);
MouseCRperc=nan(P.lastDay, size(mouseCodes,2));
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;

for m = 1 : length(mouseCodes)
    mcode = mouseCodes(m);
    ephys_dates = mice([mice.code] == mcode).ephysdates;
    ephys_dates_ = split(ephys_dates,'-');
    path_list = getBehaviorPathList(mcode);
    if ~isempty(ephys_dates)
        try
            ephys_dates_ = char(ephys_dates_(1,:,1) + ephys_dates_(1,:,2) + ephys_dates_(1,:,3));
        catch
            ephys_dates_ = ephys_dates_(1) + ephys_dates_(2) + ephys_dates_(3);
        end
        ephys_dates_trimmed = string(ephys_dates_(1,3:end,:));
        [~, folders_] = fileparts(path_list);
        path_list_ = path_list(cellfun(@(x) ~isnan(str2double(x)), folders_));
        [~, folders] = fileparts(path_list_);
        if any(contains(folders, ephys_dates_trimmed))
            path_list = path_list_(1:find(contains(folders, ephys_dates_trimmed(1)))-1);
        else
            path_list = path_list_;
        end
    end

    lastDay = min(P.lastDay, length(path_list));

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

        trial_data(f) = behavior_trial_data;
        smoothedTraces = smootheTraces(trial_data(f).eyelidpos);
        traces_new(f).traces = normalizeEyelidTraces(smoothedTraces, trial_data(f));
        trial_data(f).tracesnorm = traces_new(f).traces;
    end


    trial_data = trial_data_recalc(trial_data);

    if length(trial_data(f).CRamp)<220
        trial_data(f).CRamp(length(trial_data(f).CRamp)+1:220) = nan;
    end

    for j=1:lastDay
        Mousetraces(j,(m-1)*P.nTrials+1:m*P.nTrials,:) = trial_data(j).tracesnorm(1:P.nTrials,:);
        MouseCRamp(j,m,1:length(trial_data(j).CRamp)) = trial_data(j).CRamp';
        Mousecramp5(j,m,1:length(trial_data(j).CRamp5)) = trial_data(j).CRamp5';
        MouseCRperc(j,m) = numel(trial_data(j).CRamp5)/P.nCSplusRegTrials*100;
    end
    MouseCRperc(MouseCRperc == 0) = nan;

    type(trial_data(1).c_csdur == 0) = 0;
    type(trial_data(1).c_usdur == 0) = 2;
    Mousetype = ones(length(mouseCodes), P.nTrials);

    for i = 1 : size(Mousetype,1)
        Mousetype(i, :) = type;
    end

    ts=trial_data(1).tm(1,:);



end

for day = 1:P.lastDay
    MousetracesWithCS = Mousetraces(day, Mousetype == 1 | Mousetype == 2, :); % Extract relevant traces for the day
    MousetracesWithCS = reshape(MousetracesWithCS, [size(MousetracesWithCS, 2), size(MousetracesWithCS, 3)]);
    % Compute baseline and CR peak for all trials
    baseline = nanmean(MousetracesWithCS(:,baseMin:baseMax), 2);
    CRPeak = nanmax(MousetracesWithCS(:,CSmin:CSmax)');
    % Calculate threshold based on the amplitude from baseline to CR peak
    threshold = baseline + (CRPeak' - baseline) * P.thresCRonset / 100;
  

    for i = 1:size(MousetracesWithCS, 1)
        % Find the indices where eyelid velocity exceeds threshold
        idcs = find(MousetracesWithCS(i, CSmin:CSmax) > threshold(i), 1);
        [~, CRPeakIdx] = max(MousetracesWithCS(i, CSmin:CSmax));
        diffValues = diff(MousetracesWithCS(i, CSmin+idcs-1:CSmin+CRPeakIdx-1));
        % Check if the values are increasing after the onset index
        if ~isempty(idcs)
            if all(diffValues > 0)
                CROnsetIndices2(i) = idcs;
            else
                % Select the next index as the onset index
                nextIdx = max(find(diffValues <= 0)) + idcs;
                CROnsetIndices2(i) = nextIdx;
            end
        else
            CROnsetIndices2(i) = nan;
        end
        if CROnsetIndices2(i) < 6
            CROnsetIndices2(i) = nan;
        end
    end
    Mouseonsets(day, ~isnan(CROnsetIndices2)) = ts(CROnsetIndices2(~isnan(CROnsetIndices2)) + baseMax);
end

Mouseonsets(isnan(Mousecramp5(1,:))) = nan;

MouseCRpercstd = nanstd(MouseCRperc,1,2);
MouseCRpercSEM = MouseCRpercstd/length(mouseCodes);

for j =1:P.lastDay
    Mousecrampmean(j)=nanmean(MouseCRamp(j,:,:), 'all');
    Mousecramp5mean(j)=nanmean(Mousecramp5(j,:,:), 'all');
    Mousecrampstd(j) = nanstd(MouseCRamp(j,:,:),0, 'all');
    MousecrampSEM(j) = Mousecrampstd(j)/sqrt(length(mouseCodes));
    Mousecramp5std(j) = nanstd(Mousecramp5(j,:,:),0,'all');
    Mousemean(j,:) = nanmean(Mousetraces(j,Mousetype(1:end) == 1,:),2);
end
MouseStruct.onsets = Mouseonsets;
MouseStruct.traces = Mousetraces;
MouseStruct.cramp = MouseCRamp;
MouseStruct.cramp5 = Mousecramp5;
MouseStruct.CRperc = MouseCRperc;
MouseStruct.type = Mousetype;
MouseStruct.crampmean = Mousecrampmean;
MouseStruct.cramp5mean = Mousecramp5mean;
MouseStruct.crampstd = Mousecrampstd;
MouseStruct.crampSEM = MousecrampSEM;
MouseStruct.cramp5std = Mousecramp5std;
MouseStruct.mean = Mousemean;
MouseStruct.CRpercstd = MouseCRpercstd;
MouseStruct.CRpercSEM = MouseCRpercSEM;
MouseStruct.ts = ts;
MouseStruct.onsets = Mouseonsets;
if nargin < 2
    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s.mat", mcode)), "MouseStruct");
else
    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)), "MouseStruct");
end
end

function trial_data = trial_data_recalc(trial_data)
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
for s = 1:length(trial_data)
    traces_norm = trial_data(s).tracesnorm;
    trial_data(s).meanAll = nanmean(traces_norm);
    trial_data(s).meanCS_US = nanmean(traces_norm(trial_data(s).c_usdur~=0,:));
    trial_data(s).meanCSonly = nanmean(traces_norm(trial_data(s).c_usdur==0,:));
    trial_data(s).meanUSonly = nanmean(traces_norm(trial_data(s).c_csdur==0,:));

    %% CR amp calculation
    baseline_min(s) = nanmin(nanmean(traces_norm(:,baseMin:baseMax)));
    meanUR(s) = nanmean(nanmax(traces_norm(trial_data(s).c_usdur~=0,USmin:USmax)'));
    fullBlinkRange = meanUR(s) - baseline_min(s);

    trial_data(s).CRamp= nanmax(traces_norm(trial_data(s).c_csdur==P.CSdur,CSmin:CSmax)');
    trial_data(s).CRamp = (trial_data(s).CRamp-baseline_min(s)) / fullBlinkRange * 100;
    trial_data(s).CRamp = trial_data(s).CRamp(1:min(P.nCSplusRegTrials,length(trial_data(s).CRamp)));
    trial_data(s).CRamp5 = trial_data(s).CRamp(trial_data(s).CRamp > P.thresCRperc);

    trial_data(s).meanCR = nanmean(trial_data(s).CRamp);
    trial_data(s).stdCR = nanstd(trial_data(s).CRamp);
    trial_data(s).meanCR5 = nanmean(trial_data(s).CRamp5);
    trial_data(s).stdCR5 = nanstd(trial_data(s).CRamp5);
end
end


