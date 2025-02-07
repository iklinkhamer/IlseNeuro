%% I.K. 16-05-2024
function ephysOverviewIK()
clear
close all
app = getUImain(scriptname=mfilename);
while true
    disp("Waiting for app response...")
    uiwait(app.UIFigure)
    disp("Calculating...")
    if ~isvalid(app)
        % If the UI figure has been closed, exit the loop
        break;
    end
    appData = evalin('base', 'appData');

    if appData.enumeration == "All mice of one gene" && appData.gene == "Shank2"
        disp("Not enough Shank2 L7 mice to do this. Please pick 'All mice of both types of a gene + location' with 'Shank2' and 'KO' and try again.")
        continue
    end

    mouseCodes = appData.mousecodes;
    mouseTypes = appData.mousetypes;

    mouseCodesDefaultMice = arrayfun ...          
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
    ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
    ephysMiceMask(end-1:end) = true;
    mouseCodeEphysIdcs = contains(mouseCodes, mouseCodesDefaultMice(ephysMiceMask));

    mouseCodes = mouseCodes(mouseCodeEphysIdcs);
    mouseTypes = mouseTypes(mouseCodeEphysIdcs);

    non_valid_idcs = find(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT");
    mouseCodes(non_valid_idcs) = [];
    mouseTypes(non_valid_idcs) = [];
    if isempty(mouseCodes)
        disp("No eyelid data found. Pick another mouse.")
        continue
    end

    if appData.enumeration == "All mice of one gene" || appData.enumeration == "All mice of both types of a gene + location"
        dataStruct = struct;
        selectedMiceTypes = sort(unique(appData.mousetypes));
        types = IkUtils.getParams().types;

        for t = 1:length(selectedMiceTypes)
            t_parts = split(selectedMiceTypes(t), ' ');
            mouseTypeCC = strjoin(t_parts, '');
            idx = find([types.full] == mouseTypeCC);
            Type = types.type(idx);
            Gene = types.gene(idx);
            Loc = types.loc(idx);
            dataStruct(t).type = selectedMiceTypes(t);
            mouseListType = mouseCodes(mouseTypes == selectedMiceTypes(t));
            dataStruct(t).data = getBehaviorEyeDataEphys(mouseListType, Gene, Loc, Type);
        end

        %% Final figures
        dataStructTypes = [dataStruct.type];
        if appData.enumeration == "All mice of both types of a gene + location"
            WT = dataStruct(contains(dataStructTypes, "WT")).data;
            MUT= dataStruct(contains(dataStructTypes, "MUT")).data;
            mouseListWT = dataStruct(contains([dataStruct.type], "WT")).data.codes;
            mouseListMUT = dataStruct(contains([dataStruct.type], "MUT")).data.codes;
            makePlots(WT, mouseListWT, Gene, Loc, "WT")
            makePlots(MUT, mouseListMUT, Gene, Loc, "MUT")

        else
            WTKO = dataStruct(contains(dataStructTypes, "WT") & contains(dataStructTypes, "KO")).data;
            MUTKO = dataStruct(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "KO")).data;
            WTL7 = dataStruct(contains(dataStructTypes, "WT") & contains(dataStructTypes, "L7")).data;
            MUTL7 = dataStruct(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "L7")).data;
            mouseListWTKO = mouseListType(contains(dataStructTypes, "WT") & contains(dataStructTypes, "KO"));
            mouseListMUTKO = mouseListType(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "KO"));
            makePlots(WTKO, mouseListWTKO, "WT")
            makePlots(MUTKO, mouseListMUTKO, "MUT")

            mouseListWTL7 = mouseListType(contains(dataStructTypes, "WT") & contains(dataStructTypes, "L7"));
            mouseListMUTL7 = mouseListType(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "L7"));
            makePlots(WTL7, mouseListWTL7, "WT")
            makePlots(MUTL7, mouseListMUTL7, "MUT")
        end
    end
    disp("End of loop. Waiting for new UI input.")
end
end

function makePlots(data_eye, mouseList, Gene, Loc, Type)
p = IkUtils.getParams();
xlim_min = -0.1;
xlim_max = 0.5;
allMiceMeanTraces = [];
allMiceEphysTraces = [];
allMiceEphysFreqTrialsFac = [];
allMiceEphysFreqTrialsSup = [];
allMiceEyeTracesTrials = [];
allMiceEphysFreqTrials = [];

for m = 1:numel(mouseList)
    mcode = mouseList(m);
    miceInfo = defaultMice();
    mouseNames = [miceInfo.name];
    mname = mouseNames([miceInfo.code] == mcode);

    data = getData(mcode);
    for s = 1: length(data)
        if isempty(data(s).neuron)
            continue
        end
        modulationMasks = loadSimpleChannelModulationMasks(mcode);
        simpleMasks = loadSimpleMasks(mcode);
        modMasksSes = modulationMasks{s};
        modMasksSes = modMasksSes(simpleMasks{s});
        cs_fac_mask = [];
        cs_sup_mask = [];
        for ms = 1:length(modMasksSes)
            cs_fac_mask = [cs_fac_mask, modMasksSes{ms}.cs_facilitation];
            cs_sup_mask = [cs_sup_mask, modMasksSes{ms}.cs_suppression];
        end

        sessionDataNeurons = data(s).neuron(simpleMasks{s});
        CSfacilitatingNeuronsSession = data(s).neuron(logical(cs_fac_mask));
        CSsuppressingNeuronsSession = data(s).neuron(logical(cs_sup_mask));

        trialFreqFac = [];
        trialFreqSup = [];
        trialFreqAll = [];

        for n2 = 1:length(CSsuppressingNeuronsSession)
            rasterTrialsSup = CSsuppressingNeuronsSession(n2).RasterXY_cs_filtered(2,1:3:end);
            uniqueTrials = unique(sort(rasterTrialsSup));
            trialFreq = [];
            for t = 1:220
                if ~ismember(t, uniqueTrials)
                    trialFreq = [trialFreq;nan];
                else
                    spikes = CSsuppressingNeuronsSession(n2).RasterXY_cs_filtered(1,1:3:end);
                    trialspikes = spikes(rasterTrialsSup == t);
                    range = p.sspkRanges.cs.min:p.BinW:p.sspkRanges.cs.max;
                    trialspikecounts = histcounts(trialspikes, range);
                    trialspikerate = sum(trialspikecounts)/(p.sspkRanges.cs.max-p.sspkRanges.cs.min);
                    trialFreq = [trialFreq;trialspikerate];
                end
            end
            trialFreqSup = [trialFreqSup, {trialFreq}];
        end
        for n2 = 1:length(CSfacilitatingNeuronsSession)
            rasterTrialsFac = CSfacilitatingNeuronsSession(n2).RasterXY_cs_filtered(2,1:3:end);
            uniqueTrials = unique(sort(rasterTrialsFac));
            trialFreq = [];
            for t = 1:220
                if ~ismember(t, uniqueTrials)
                    trialFreq = [trialFreq;nan];
                else                    
                    spikes = CSfacilitatingNeuronsSession(n2).RasterXY_cs_filtered(1,1:3:end);
                    trialspikes = spikes(rasterTrialsFac == t);
                    trialspikecounts = histcounts(trialspikes, p.sspkRanges.cs.min:p.BinW:p.sspkRanges.cs.max);
                    trialspikerate = sum(trialspikecounts)/(p.sspkRanges.cs.max-p.sspkRanges.cs.min);
                    trialFreq = [trialFreq;trialspikerate];
                end
            end
            trialFreqFac = [trialFreqFac, {trialFreq}];
        end
        for n2 = 1:length(sessionDataNeurons)
            rasterTrials = sessionDataNeurons(n2).RasterXY_cs_filtered(2,1:3:end);
            uniqueTrials = unique(sort(rasterTrials));
            trialFreq = [];
            for t = 1:220
                if ~ismember(t, uniqueTrials)
                    trialFreq = [trialFreq;nan];
                else                    
                    spikes = sessionDataNeurons(n2).RasterXY_cs_filtered(1,1:3:end);
                    trialspikes = spikes(rasterTrials == t);
                    trialspikecounts = histcounts(trialspikes, p.sspkRanges.cs.min:p.BinW:p.sspkRanges.cs.max);
                    trialspikerate = sum(trialspikecounts)/(p.sspkRanges.cs.max-p.sspkRanges.cs.min);
                    trialFreq = [trialFreq;trialspikerate];
                end
            end
            trialFreqAll = [trialFreqAll, {trialFreq}];
        end

        if isempty(sessionDataNeurons)
            continue
        end

        trial_type_mask_csus_ = data_eye.type(:,1:p.nTrials) == 1;
        trial_type_mask_csus = [];
        for ttm = 1:size(trial_type_mask_csus_, 1)
            trial_type_mask_csus = [trial_type_mask_csus, trial_type_mask_csus_(ttm,:)];
        end
        sessionMask = false(1, size(data_eye.traces,2));
        sessionMask((m-1)*p.nTrials+1:m*p.nTrials) = true;
        mask = trial_type_mask_csus & sessionMask;
        traces = squeeze(data_eye.traces(s, mask, :));
        meantraces = nanmean(traces);

        if all(isnan(meantraces))
            continue
        end

        ephysTracesSession = [];
        for n = 1:length(sessionDataNeurons)
            ephysTracesSession = [ephysTracesSession; sessionDataNeurons(n).psth_cs_reset_filtered(2,:)];
        end

            trial_type_mask_csus_csonly_ = data_eye.type(:,1:p.nTrials) == 1 | data_eye.type(:,1:p.nTrials) == 0;
        trial_type_mask_csus_csonly = [];
        for ttm = 1:size(trial_type_mask_csus_csonly_, 1)
            trial_type_mask_csus_csonly = [trial_type_mask_csus_csonly, trial_type_mask_csus_csonly_(ttm,:)];
        end
        mask_csus_csonly = trial_type_mask_csus_csonly & sessionMask;
        traces_ephys_trials = squeeze(data_eye.traces(s, mask_csus_csonly, :));
        meantraces_ephys_trials = nanmean(traces_ephys_trials);

        allMiceEphysTraces = [allMiceEphysTraces; {ephysTracesSession}];
        allMiceMeanTraces = [allMiceMeanTraces; meantraces_ephys_trials];
        allMiceEphysFreqTrialsFac = [allMiceEphysFreqTrialsFac; {trialFreqFac}];
        allMiceEphysFreqTrialsSup = [allMiceEphysFreqTrialsSup; {trialFreqSup}];
        allMiceEyeTracesTrials = [allMiceEyeTracesTrials; {traces_ephys_trials}];
        allMiceEphysFreqTrials = [allMiceEphysFreqTrials; {trialFreqAll}];

        if true
            continue
        end

        figure('Color', 'white','Position', [10 0 1000 1000]);     % Create a new figure for each session
        subplot('Position', [0.1, 0.45, 0.8, 0.35]); % Plot the mean trace in the upper subplot
        hold on;
        % plot blue and green blocks
        xlim([xlim_min, xlim_max]);
        set(gca, 'tickdir', 'out');

        x1 = 0;
        x2 = p.sspkRanges.cs.max;
        y1 = min(meantraces);
        y2 = max(max(traces));
        patch ...
            ( [x1 x1 x2 x2] ...
            , [y1 y2 y2 y1] ...
            , getColors().histogram.csRegion ...
            , EdgeColor = 'none' ...
            , FaceAlpha = 0.1 ...
            )

        x1 = p.sspkRanges.us.min;
        x2 = x1 + p.USdur/1000;
        patch ...
            ( [x1 x1 x2 x2] ...
            , [y1 y2 y2 y1] ...
            , getColors().histogram.usRegion ...
            , EdgeColor = 'none' ...
            , FaceAlpha = 0.1 ...
            )
        plot(data_eye.ts(1:41)/1000, traces(:,1:41), 'Color', [0.5, 0.5, 0.5]);
        plot(data_eye.ts(41:91)/1000, traces(:,41:91), 'Color', [0 0 1]);
        plot(data_eye.ts(91:95)/1000, traces(:,91:95), 'Color', [0, 0.5, 0]);
        plot(data_eye.ts(95:p.nTimeSteps)/1000, traces(:,95:p.nTimeSteps), 'Color', [0.5, 0.5, 0.5]);
        plot(data_eye.ts/1000, meantraces, 'LineWidth', 3, 'Color', [0,0,0]);
        title(sprintf('%s %s - Session %d', mname, mcode, s));
        %         xticklabels({});
        xticks([]);
        %         xticklabels([0, 250]);
        ylim([y1, y2])
        yticks([y1 y2])
        yticklabels([0 1])
        ylabel("Eyelid Closure")


        arrayfun ...
            ( @(eventTime) ...
            plot ...
            ( [eventTime eventTime], [y1, y2] ...
            , '--' ...
            , color = [.5 .5 .5] ...
            ) ...
            , p.delayEventsCs ...
            );

        box('off')

        %% Raster
        maxHeightOffset = 0;
        subplot('Position', [0.1, 0.1, 0.8, 0.35]);
        hold on


        for n =1:length(sessionDataNeurons)
            if isempty(sessionDataNeurons(n).RasterXY_cs_filtered)
                continue
            end
            csIndices = sessionDataNeurons(n).RasterXY_cs_filtered(3, :) == 0;
            csIndices(find(csIndices == 1)+1) = true;
            usIndices = sessionDataNeurons(n).RasterXY_cs_filtered(3, :) == 2;
            usIndices(find(usIndices == 1) + 1) = true;

            rasterXY = sessionDataNeurons(n).RasterXY_cs_filtered;
            rasterXY_filtered = rasterXY;
            rasterXY_filtered(:, csIndices|usIndices) = [];

            uniqueTrials = unique(rasterXY_filtered(2, 1:3:end));  % Unique trial values
            adapted_height = 1:numel(uniqueTrials);      % Adapted heights for each unique trial

            trialIndices = ismember(rasterXY_filtered(2, :), uniqueTrials); % Find indices of occurrences of unique trials in rasterXY

            % % Find indices of occurrences of unique trials in rasterXY
            [~, ~, trialIdx] = unique(rasterXY_filtered(2, 1:3:end));  % Extract unique trial indices

            rasterXY_adapted = rasterXY_filtered; % Assign adapted heights directly
            rasterXY_adapted(2, trialIndices) = adapted_height(trialIdx);

            nextIndices = false(size(trialIndices)); % Create a logical index shifted by one position
            nextIndices(1:end-2) = trialIndices(3:end);
            nextIndices(end-1) = true;

            rasterXY_adapted(2, nextIndices) = rasterXY_adapted(2, trialIndices) + p.raster_spike_height; % Increment the following indices by raster_spike_height

            plot ...
                ( rasterXY_adapted(1,:) * 1000 ...
                , rasterXY_adapted(2,:) + maxHeightOffset...
                , Color = [.2 .2 .2] ...
                );

            maxHeightOffset = floor(max(rasterXY_adapted(2,:))) + maxHeightOffset + 1;
        end
        eventMn = 0;

        arrayfun ...
            ( @(eventTime) ...
            plot ...
            ( [eventTime eventTime] * 1000, [eventMn, maxHeightOffset] ...
            , '--' ...
            , color = [.5 .5 .5] ...
            ) ...
            , p.delayEventsCs ...
            );

        ylim([0 maxHeightOffset]);
        xlim([xlim_min*1000 xlim_max*1000])
        xticks([0 250]);
        set(gca, 'tickdir', 'out');
        xticklabels([0 250]);
        xlabel("Time (ms)")

        yticks([])
        ylabel("Trials all neurons")

        box('off')

        % plot blue and green blocks
        x1 = 0;
        x2 = p.sspkRanges.cs.max * 1000;
        y1 = eventMn;
        y2 = y1 + maxHeightOffset;
        patch ...
            ( [x1 x1 x2 x2] ...
            , [y1 y2 y2 y1] ...
            , getColors().histogram.csRegion ...
            , EdgeColor = 'none' ...
            , FaceAlpha = 0.1 ...
            )

        x1 = p.sspkRanges.us.min * 1000;
        x2 = p.sspkRanges.us.min* 1000 + p.USdur;

        patch ...
            ( [x1 x1 x2 x2] ...
            , [y1 y2 y2 y1] ...
            , getColors().histogram.usRegion ...
            , EdgeColor = 'none' ...
            , FaceAlpha = 0.1 ...
            )

        fname = sprintf('ephys-beh-overview_%s_S%d.png', mcode, s);
        file = fullfile(p.figPath, "Ephys and Behavior Overview",fname);
        saveas(gca,file);
        fname = sprintf('ephys-beh-overview_%s_S%d.eps', mcode, s);
        file = fullfile(p.figPath, "Ephys and Behavior Overview",fname);
        print(gcf, '-depsc', '-painters', file)


        %         pause();
    end
end
correlationData = struct;
correlationData.eye_traces = allMiceEyeTracesTrials;
correlationData.mean_eye_traces = allMiceMeanTraces;
correlationData.mean_spike_rate = allMiceEphysTraces;
correlationData.freq_fac_trials = allMiceEphysFreqTrialsFac;
correlationData.freq_sup_trials = allMiceEphysFreqTrialsSup;
correlationData.freq_trials = allMiceEphysFreqTrials;

save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("correlationData_%s%s%s.mat", Gene, Loc, Type)), "correlationData");

% load('/home/i.klinkhamer/Documents/Data/other/ephysOverviewData_Shank2KOWT.mat')
% load('/home/i.klinkhamer/Documents/Data/other/correlationData_Shank2KOWT.mat')
p = IkUtils.getParams();
eyetraceCSrange = find(data_eye.ts >=0 & data_eye.ts < 250);
bins = -p.t_raster_edges : p.spikeTimesBinW : p.t_raster_edges;
range = find(bins >= 0 & bins < p.sspkRanges.cs.max);
baseMin = p.tracesRanges.baseline.min;
baseMax = p.tracesRanges.baseline.max;
USmin = p.tracesRanges.us.min;
USmax = p.tracesRanges.us.max;

group = Gene+Loc+Type;

mean_spike_rate_neuron = [];
eye_trace_all = [];
for cd = 1:length(correlationData.freq_trials)
    mean_eye_trace(cd,:) = correlationData(1).mean_eye_traces(cd,:);
    mean_eye_trace_cs(cd, :) = mean_eye_trace(cd, eyetraceCSrange);
    % Define the factor by which you want to upsample
    upsample_factor = 10;
    % Create the original time points
    original_time = linspace(1, size(mean_eye_trace_cs(cd,:),2), size(mean_eye_trace_cs(cd,:),2));
    % Create the new time points after upsampling
    new_time = linspace(1, size(mean_eye_trace_cs(cd,:),2), size(mean_eye_trace_cs(cd,:),2)*upsample_factor);
    % Interpolate the data
    mean_eye_trace_cs_upsampled(cd,:) = interp1(original_time, mean_eye_trace_cs(cd,:).', new_time, 'spline');

    mean_spike_rate(cd) = [correlationData(1).mean_spike_rate(cd,:)];
    mean_spike_rate_neuron_ = mean_spike_rate{cd};
    for n3 = 1:size(mean_spike_rate_neuron_,1)
        eye_trace_all = [eye_trace_all; mean_eye_trace_cs_upsampled(cd,:)];
    end
    mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_];
end
mean_spike_rate_cs = mean_spike_rate_neuron(:,range);

% Compute normalized cross-correlation
for cc = 1:size(eye_trace_all,1)
    [c(cc,:), lags(cc,:)] = xcorr(mean_spike_rate_cs(cc,:), eye_trace_all(cc, :), 'coeff');
end

% % Plot the normalized cross-correlation
% figure;
% stem(lags, c);
% xlabel('Lag');
% ylabel('Normalized Cross-correlation');
% ylim([0 1])
% title('Normalized Cross-correlation between x and y');

% Plot using plot function instead of stem
figure; hold on 
for cc2 = 1:size(eye_trace_all, 1)
plot(lags(cc2,:), c(cc2,:)); 
end
x = [-500:500];
plot(x, zeros(size(x,1), size(x,2)), '--', 'Color', 'k');
plot(zeros(size(x,1), size(x,2)), x, '--', 'Color', 'k');
xlabel('Lag');
ylabel('Cross-correlation');
ylim([-1 1])
title('Cross-correlation between average spike rate and eye trace');
hold off

fname = sprintf('xcorr_spikes_eye_CS_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('xcorr_spikes_eye_CS_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)

fac_mask = [];
freq_fac_trials = [correlationData.freq_fac_trials];
for fft = 1: length(correlationData.freq_fac_trials)
if ~isempty(freq_fac_trials{fft})
    fac_mask = [fac_mask; 1];
else
    fac_mask = [fac_mask; 0];
end
end
fac_mask = logical(fac_mask);
fac_idcs = find(fac_mask);

sup_mask = [];
freq_sup_trials = [correlationData.freq_sup_trials];
for fft = 1: length(correlationData.freq_sup_trials)
if ~isempty(freq_sup_trials{fft})
    sup_mask = [sup_mask; 1];
else
    sup_mask = [sup_mask; 0];
end
end
sup_mask = logical(sup_mask);
sup_idcs = find(sup_mask);

freq_neurons_all = [];
CRamp_CS_all = [];
for fn = 1:length(fac_idcs)
    freq_neuron_ = freq_fac_trials{fac_idcs(fn)};
    for n = 1:size(freq_neuron_,2)
        freq_neuron = freq_neuron_{n};
        freq_neurons_all = [freq_neurons_all, freq_neuron'];
    end

    eye_traces_neuron = correlationData.eye_traces{fac_idcs(fn)};
    CRamp_CS = max(eye_traces_neuron(:, eyetraceCSrange)');
    baseline_min = nanmin(nanmean(eye_traces_neuron(:,baseMin:baseMax)));
    meanUR = nanmean(nanmax(eye_traces_neuron(:,USmin:USmax)'));
    fullBlinkRange = meanUR - baseline_min;

    CRamp_CS = (CRamp_CS-baseline_min) / fullBlinkRange * 100;

    CRamp_CS(CRamp_CS < p.thresCRperc) = nan;
    CR_temp_var = CRamp_CS;
    for n = 1:(size(freq_neuron_,2)-1)
        CRamp_CS = [CRamp_CS, CR_temp_var];
    end
    CRamp_CS_all = [CRamp_CS_all, CRamp_CS];

end


mask = ~isnan(freq_neurons_all) & ~isnan(CRamp_CS_all);

% Calculate the correlation coefficient
[R, P] = corrcoef(freq_neurons_all(mask), CRamp_CS_all(mask));
correlation_coefficient = R(1, 2); % Extract the correlation coefficient

% Plot the scatter plot
figure; hold on
scatter(freq_neurons_all(mask), CRamp_CS_all(mask), 'filled'); % Use 'filled' to fill the markers
xlabel('frequency');
ylabel('CR amplitude');
title('Scatter Plot of rate and CR amplitude CS facilitation');
% Perform linear regression
coefficients = polyfit(freq_neurons_all(mask), CRamp_CS_all(mask), 1); % 1 indicates a linear fit

% Evaluate the fit
% x_fit = linspace(min(freq_neurons_all(mask)), max(freq_neurons_all(mask)), 1000); % Generate x values for the fit line
y_fit = polyval(coefficients, freq_neurons_all(mask)); % Compute the corresponding y values
plot(freq_neurons_all(mask), y_fit, '-r', 'LineWidth', 2);

% Display the correlation coefficient on the plot
text(min(freq_neurons_all(mask)), max(CRamp_CS_all(mask)), ['Correlation: ', num2str(correlation_coefficient), '   P: ', num2str(P(1,2))], 'FontSize', 12, 'Color', 'red');
% text(min(freq_neurons_all(mask))+100, max(CRamp_CS_all(mask)), ['P: ', num2str(P(1,2))], 'FontSize', 12, 'Color', 'red');
hold off

fname = sprintf('corr_meanratetrials_maxamp_facilitating_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('corr_meanratetrials_maxamp_facilitating_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)


%% suppressing
if ~isempty(sup_idcs)
    freq_neurons_all_sup = [];
    CRamp_CS_all_sup = [];
    for fn = 1:length(sup_idcs)
        freq_neuron_ = freq_sup_trials{sup_idcs(fn)};
        for n = 1:size(freq_neuron_,2)
            freq_neuron = freq_neuron_{n};
            freq_neurons_all_sup = [freq_neurons_all_sup, freq_neuron'];
        end

        eye_traces_neuron = correlationData.eye_traces{sup_idcs(fn)};
        CRamp_CS = max(eye_traces_neuron(:, eyetraceCSrange)');
        baseline_min = nanmin(nanmean(eye_traces_neuron(:,baseMin:baseMax)));
        meanUR = nanmean(nanmax(eye_traces_neuron(:,USmin:USmax)'));
        fullBlinkRange = meanUR - baseline_min;

        CRamp_CS = (CRamp_CS-baseline_min) / fullBlinkRange * 100;

        CRamp_CS(CRamp_CS < p.thresCRperc) = nan;
        CR_temp_var = CRamp_CS;
        for n = 1:(size(freq_neuron_,2)-1)
            CRamp_CS = [CRamp_CS, CR_temp_var];
        end
        CRamp_CS_all_sup = [CRamp_CS_all_sup, CRamp_CS];

    end

    mask = ~isnan(freq_neurons_all_sup) & ~isnan(CRamp_CS_all_sup);

    % Calculate the correlation coefficient
    [R, P] = corrcoef(freq_neurons_all_sup(mask), CRamp_CS_all_sup(mask));
    correlation_coefficient = R(1, 2); % Extract the correlation coefficient

    % Plot the scatter plot
    figure; hold on
    scatter(freq_neurons_all_sup(mask), CRamp_CS_all_sup(mask), 'filled'); % Use 'filled' to fill the markers
    xlabel('frequency');
    ylabel('CR amplitude');
    title('Scatter Plot of rate and CR amplitude CS suppression');
%     % Perform linear regression
%     coefficients = polyfit(freq_neurons_all_sup(mask), CRamp_CS_all_sup(mask), 1); % 1 indicates a linear fit
% 
%     % Evaluate the fit
%     x_fit = linspace(min(freq_neurons_all_sup(mask)), max(freq_neurons_all_sup(mask)), 1000); % Generate x values for the fit line
%     y_fit = polyval(coefficients, x_fit); % Compute the corresponding y values
%     plot(x_fit, y_fit, '-r', 'LineWidth', 2);

    % Perform linear regression
    coefficients = polyfit(freq_neurons_all_sup(mask), CRamp_CS_all_sup(mask), 1); % 1 indicates a linear fit

    % Evaluate the fit
%     x_fit = linspace(min(freq_neurons_all_sup(mask)), max(freq_neurons_all_sup(mask)), 1000); % Generate x values for the fit line
    y_fit = polyval(coefficients, freq_neurons_all_sup(mask)); % Compute the corresponding y values
    plot(freq_neurons_all_sup(mask), y_fit, '-r', 'LineWidth', 2);

    % Display the correlation coefficient on the plot
    text(min(freq_neurons_all_sup(mask)), max(CRamp_CS_all_sup(mask)), ['Correlation: ', num2str(correlation_coefficient), '   P: ', num2str(P(1,2))], 'FontSize', 12, 'Color', 'red');
    % text(min(freq_neurons_all(mask))+100, max(CRamp_CS_all(mask)), ['P: ', num2str(P(1,2))], 'FontSize', 12, 'Color', 'red');
    hold off

    fname = sprintf('corr_meanratetrials_maxamp_suppressing_%s.png', group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    saveas(gca,file);
    fname = sprintf('corr_meanratetrials_maxamp_suppressing_%s.eps', group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    print(gcf, '-depsc', '-painters', file)
end

end
