%% I.K. 16-05-2024
function ephysOverview2IK()
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
            [optimal_lags_fac_WT, optimal_lags_sup_WT] = makePlots(WT, mouseListWT, Gene, Loc, "WT");
            [optimal_lags_fac_MUT, optimal_lags_sup_MUT] = makePlots(MUT, mouseListMUT, Gene, Loc, "MUT");

            [p_fac, h_fac, stats_fac] = ranksum(optimal_lags_fac_WT, optimal_lags_fac_MUT);
            disp(['p-value facilitating: ', num2str(p_fac), ' ', Gene, Loc]);
            try
                [p_sup, h_sup, stats_sup] = ranksum(optimal_lags_sup_WT, optimal_lags_sup_MUT);
                disp(['p-value suppressing: ', num2str(p_sup), ' ', Gene, Loc]);
            catch
            end



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

function [optimal_lags_no_outliers_fac, optimal_lags_no_outliers_sup] = makePlots(data_eye, mouseList, Gene, Loc, Type)
skipOverviewPlots = true;
skipCorrelationPlots = false;
redoCorrelationData = false;

p = IkUtils.getParams();
group = Gene+Loc+Type;

xlim_min = -0.1;
xlim_max = 0.5;
allMiceMeanEyeTracesCS = [];

allMiceMeanEphysTraces = [];
allMiceEphysFreqTrialsFacCS = [];
allMiceEphysFreqTrialsSupCS = [];

allMiceEyeTracesTrials = [];
allMiceEphysFreqTrialsCS = [];

allMiceEphysFreqTrialsUS = [];
allMiceEphysFreqTrialsFacUS = [];
allMiceEphysFreqTrialsSupUS = [];
allMiceMeanEyeTracesUS = [];

allMiceModulationMasks = [];
allMiceCodes = [];
allMiceSessions = [];
for m = 1:numel(mouseList)
    if ~redoCorrelationData && skipOverviewPlots
        continue
    end
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
        sessionDataNeurons = data(s).neuron(simpleMasks{s});
        if isempty(sessionDataNeurons)
            continue
        end
        allMiceCodes = [allMiceCodes; mcode];
        allMiceSessions = [allMiceSessions; s];
        % Eye traces
        sessionMask = false(1, size(data_eye.traces,2));
        sessionMask((m-1)*p.nTrials+1:m*p.nTrials) = true;


        trial_type_mask_csus_ = data_eye.type(:,1:p.nTrials) == 1;
        trial_type_mask_csus = [];
        for ttm = 1:size(trial_type_mask_csus_, 1)
            trial_type_mask_csus = [trial_type_mask_csus, trial_type_mask_csus_(ttm,:)];
        end
        mask_csus = trial_type_mask_csus & sessionMask;
        traces_ephys_trials_CSUS_combined = squeeze(data_eye.traces(s, mask_csus, :));
        meantraces_ephys_trials_CSUS_combined = nanmean(traces_ephys_trials_CSUS_combined);

        if all(isnan(meantraces_ephys_trials_CSUS_combined))
            continue
        end

        if redoCorrelationData || ~exist(sprintf("/home/i.klinkhamer/Documents/Data/other/correlationData_%s.mat", group), 'file')
            modMasksSes = modulationMasks{s};
            modMasksSes = modMasksSes(simpleMasks{s});
            cs_fac_mask = [];
            cs_sup_mask = [];
            us_fac_mask = [];
            us_sup_mask = [];
            for ms = 1:length(modMasksSes)
                cs_fac_mask = [cs_fac_mask, modMasksSes{ms}.cs_facilitation];
                cs_sup_mask = [cs_sup_mask, modMasksSes{ms}.cs_suppression];
                us_fac_mask = [us_fac_mask, modMasksSes{ms}.us_facilitation];
                us_sup_mask = [us_sup_mask, modMasksSes{ms}.us_suppression];
            end

            CSfacilitatingNeuronsSession = data(s).neuron(logical(cs_fac_mask));
            CSsuppressingNeuronsSession = data(s).neuron(logical(cs_sup_mask));
            USfacilitatingNeuronsSession = data(s).neuron(logical(us_fac_mask));
            USsuppressingNeuronsSession = data(s).neuron(logical(us_sup_mask));

            trialFreqFacCS = calcAveSpkRateROI(CSfacilitatingNeuronsSession, p.sspkRanges.cs);
            trialFreqSupCS = calcAveSpkRateROI(CSsuppressingNeuronsSession, p.sspkRanges.cs);
            trialFreqAllCS = calcAveSpkRateROI(sessionDataNeurons, p.sspkRanges.cs);

            trialFreqFacUS = calcAveSpkRateROI(USfacilitatingNeuronsSession, p.sspkRanges.us);
            trialFreqSupUS = calcAveSpkRateROI(USsuppressingNeuronsSession, p.sspkRanges.us);
            trialFreqAllUS = calcAveSpkRateROI(sessionDataNeurons, p.sspkRanges.us);


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
            traces_ephys_trials_CS = squeeze(data_eye.traces(s, mask_csus_csonly, :));
            meantraces_ephys_trials_CS = nanmean(traces_ephys_trials_CS);

            allMiceMeanEphysTraces = [allMiceMeanEphysTraces; {ephysTracesSession}];
            allMiceMeanEyeTracesCS = [allMiceMeanEyeTracesCS; meantraces_ephys_trials_CS];
            allMiceMeanEyeTracesUS = [allMiceMeanEyeTracesUS; meantraces_ephys_trials_CSUS_combined];
            allMiceEphysFreqTrialsFacCS = [allMiceEphysFreqTrialsFacCS; {trialFreqFacCS}];
            allMiceEphysFreqTrialsSupCS = [allMiceEphysFreqTrialsSupCS; {trialFreqSupCS}];
            allMiceEphysFreqTrialsFacUS = [allMiceEphysFreqTrialsFacUS; {trialFreqFacUS}];
            allMiceEphysFreqTrialsSupUS = [allMiceEphysFreqTrialsSupUS; {trialFreqSupUS}];
            allMiceEyeTracesTrials = [allMiceEyeTracesTrials; {traces_ephys_trials_CS}];
            allMiceEphysFreqTrialsCS = [allMiceEphysFreqTrialsCS; {trialFreqAllCS}];
            allMiceEphysFreqTrialsUS = [allMiceEphysFreqTrialsUS; {trialFreqAllUS}];
            allMiceModulationMasks = [allMiceModulationMasks; {modMasksSes}];

        end

        if skipOverviewPlots
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
        y1 = min(meantraces_ephys_trials_CSUS_combined);
        y2 = max(max(traces_ephys_trials_CSUS_combined));
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
        plot(data_eye.ts(1:41)/1000, traces_ephys_trials_CSUS_combined(:,1:41), 'Color', [0.5, 0.5, 0.5]);
        plot(data_eye.ts(41:91)/1000, traces_ephys_trials_CSUS_combined(:,41:91), 'Color', [0 0 1]);
        plot(data_eye.ts(91:95)/1000, traces_ephys_trials_CSUS_combined(:,91:95), 'Color', [0, 0.5, 0]);
        plot(data_eye.ts(95:p.nTimeSteps)/1000, traces_ephys_trials_CSUS_combined(:,95:p.nTimeSteps), 'Color', [0.5, 0.5, 0.5]);
        plot(data_eye.ts/1000, meantraces_ephys_trials_CSUS_combined, 'LineWidth', 3, 'Color', [0,0,0]);
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
            if n~=length(sessionDataNeurons)
                plot([xlim_min*1000 xlim_max*1000], [floor(max(rasterXY_adapted(2,:)))+ maxHeightOffset, floor(max(rasterXY_adapted(2,:)))+ maxHeightOffset], '-', 'Color', 'k')
            end
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

if skipCorrelationPlots
    return
end

if redoCorrelationData || ~exist(sprintf("/home/i.klinkhamer/Documents/Data/other/correlationData_%s.mat", group), 'file')
    correlationData = struct;
    correlationData.eye_traces = allMiceEyeTracesTrials;
    correlationData.mean_spike_rate = allMiceMeanEphysTraces;
    correlationData.mod_masks = allMiceModulationMasks;
    correlationData.m_codes = allMiceCodes;
    correlationData.sessions = allMiceSessions;

    correlationData.cs.mean_eye_traces = allMiceMeanEyeTracesCS;
    correlationData.cs.freq_fac_trials = allMiceEphysFreqTrialsFacCS;
    correlationData.cs.freq_sup_trials = allMiceEphysFreqTrialsSupCS;
    correlationData.cs.freq_trials = allMiceEphysFreqTrialsCS;

    correlationData.us.mean_eye_traces = allMiceMeanEyeTracesUS;
    correlationData.us.freq_fac_trials = allMiceEphysFreqTrialsFacUS;
    correlationData.us.freq_sup_trials = allMiceEphysFreqTrialsSupUS;
    correlationData.us.freq_trials = allMiceEphysFreqTrialsUS;

    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("correlationData_%s%s%s.mat", Gene, Loc, Type)), "correlationData");
else
    load(sprintf('/home/i.klinkhamer/Documents/Data/other/correlationData_%s.mat', group))
end
% load('/home/i.klinkhamer/Documents/Data/other/ephysOverviewData_Tsc1KOWT.mat')
% load('/home/i.klinkhamer/Documents/Data/other/correlationData_Tsc1KOWT.mat')
p = IkUtils.getParams();

bins = -p.t_raster_edges : p.spikeTimesBinW : p.t_raster_edges;

rangeCS = find(bins >= 0 & bins < p.sspkRanges.cs.max);
 rangeBase = find(bins < 0);
 test = correlationData.mean_spike_rate{1};
  testBase = test(:,rangeBase);
 meanBase = mean(testBase,2);

[eye_trace_fac_cs, mean_spike_rate_cs_fac] = calcStatsCC(correlationData, "cs", rangeCS, "fac");
rangeUS = find(bins >= p.sspkRanges.us.min & bins < p.sspkRanges.us.max);
[eye_trace_sup_cs, mean_spike_rate_cs_sup] = calcStatsCC(correlationData, "cs", rangeCS, "sup");

for cc = 1:size(eye_trace_fac_cs,1) % Compute normalized cross-correlation for facilitation
    [c(cc,:), lags(cc,:)] = xcorr(eye_trace_fac_cs(cc, :), mean_spike_rate_cs_fac(cc,:), 'coeff');
    [~, max_idx(cc)] = max(c(cc,:)); % Find the optimal lag (lag with the highest correlation)
    optimal_lags(cc) = lags(cc,max_idx(cc));        
    correlation_coeffs(cc) = c(cc,max_idx(cc)); % Store the maximum correlation coefficient
end
      
horizontal_jitter = @(x) x + (rand(size(x)) - 0.5) * 0.1; % Horizontal jittering function % Adjust the multiplier for more or less horizontal jitter

figure; % Plot the correlation coefficients using a boxplot
boxplot(optimal_lags);
outliers = get(findobj(gca, 'Tag', 'Outliers'), 'YData');
optimal_lags_no_outliers_fac = optimal_lags(~ismember(optimal_lags, outliers))/2;

[p_corr, ~, ~] = signrank(correlation_coeffs(~ismember(optimal_lags, outliers)));% Perform Wilcoxon signed-rank test on the correlation coefficients
[p_lag, ~, ~] = signrank(optimal_lags_no_outliers_fac);% Perform Wilcoxon signed-rank test on the optimal lags
sd_lags = std(optimal_lags_no_outliers_fac);
hold on;
plot(horizontal_jitter(ones(size(optimal_lags))), optimal_lags, 'r.', 'MarkerSize', 10); % Add horizontally jittered points as red dots
xlabel(sprintf('%s  %s', group, 'Facilitation'));
ylim([-500 500]);
yticks(-500:100:500);
yticklabels(-250:50:250);
ylabel('Optimal Lag (ms)');
title(['Optimal Lags (p = ', num2str(p_lag), ')']);
set(gca, 'tickdir', 'out');
box('off');
hold off;
disp (["Facilitating neurons", group]);
disp(['p-value for correlation coefficients: ', num2str(p_corr)]);
disp(['p-value for optimal lags: ', num2str(p_lag)]);

disp(['Optimal lag: ', num2str(mean(optimal_lags_no_outliers_fac)), ' ms']);
disp(['SD lags: ', num2str(sd_lags), 'ms']);

fname = sprintf('boxplot_wil_CS_facilitation_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('boxplot_wil_CS_facilitation_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)


figure; hold on
for cc2 = 1:size(eye_trace_fac_cs, 1)
    plot(lags(cc2,:), c(cc2,:));
end
plot(mean(lags), mean(c), 'LineWidth',2, 'color', 'k');
x = [-500:500];
xticks(-500:100:500)
xticklabels(-250:50:250)
plot(x, zeros(size(x,1), size(x,2)), '--', 'Color', 'k');
plot(zeros(size(x,1), size(x,2)), x, '--', 'Color', 'k');
xlabel('Lag');
ylabel('Cross-correlation');
ylim([-1 1])
title('Cross-correlation between average spike rate and eye trace during CS for facilitating neurons');
set(gca, 'tickdir', 'out');
box('off');
hold off

fname = sprintf('xcorr_spikes_eye_CS_facilitation_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('xcorr_spikes_eye_CS_facilitation_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)


c = []; lags = []; optimal_lags = []; correlation_coeffs = []; % Compute normalized cross-correlation for suppression
for cc = 1:size(eye_trace_sup_cs,1)
    [c(cc,:), lags(cc,:)] = xcorr(eye_trace_sup_cs(cc, :), mean_spike_rate_cs_sup(cc,:), 'coeff');
    [~, min_idx(cc)] = min(c(cc,:)); % Find the optimal lag (lag with the highest correlation)
    optimal_lags(cc) = lags(cc,min_idx(cc));        
    correlation_coeffs(cc) = c(cc,min_idx(cc)); % Store the maximum correlation coefficient
end
      
try

horizontal_jitter = @(x) x + (rand(size(x)) - 0.5) * 0.1; % Horizontal jittering function % Adjust the multiplier for more or less horizontal jitter
figure; % Plot the correlation coefficients using a boxplot
boxplot(optimal_lags);
outliers = get(findobj(gca, 'Tag', 'Outliers'), 'YData');
optimal_lags_no_outliers_sup = optimal_lags(~ismember(optimal_lags, outliers))/2;

[p_corr, ~, ~] = signrank(correlation_coeffs(~ismember(optimal_lags, outliers)));% Perform Wilcoxon signed-rank test on the correlation coefficients
[p_lag, ~, ~] = signrank(optimal_lags_no_outliers_sup);% Perform Wilcoxon signed-rank test on the optimal lags
sd_lags = std(optimal_lags_no_outliers_sup);

hold on;
plot(horizontal_jitter(ones(size(optimal_lags))), optimal_lags, 'r.', 'MarkerSize', 10); % Add horizontally jittered points as red dots
xlabel(sprintf('%s  %s', group, 'Suppression'));
ylim([-500 500]);
yticks(-500:100:500);
yticklabels(-250:50:250);
ylabel('Optimal Lag (ms)');
title(['Optimal Lags (p = ', num2str(p_lag), ')']);
box('off');
set(gca, 'tickdir', 'out');
hold off;

disp (["Suppressing neurons", group]);
disp(['p-value for correlation coefficients: ', num2str(p_corr)]);
disp(['p-value for optimal lags: ', num2str(p_lag)]);

disp(['Optimal lag: ', num2str(mean(optimal_lags_no_outliers_sup)), ' ms']);
disp(['SD lags: ', num2str(sd_lags), 'ms']);

fname = sprintf('boxplot_wil_CS_suppression_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('boxplot_wil_CS_suppression_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)

catch
end

figure; hold on
for cc2 = 1:size(eye_trace_sup_cs, 1)
    plot(lags(cc2,:), c(cc2,:));
end
plot(mean(lags), mean(c), 'LineWidth',2, 'color', 'k');
x = [-500:500];
xticks(-500:100:500)
xticklabels(-250:50:250)
plot(x, zeros(size(x,1), size(x,2)), '--', 'Color', 'k');
plot(zeros(size(x,1), size(x,2)), x, '--', 'Color', 'k');
xlabel('Lag');
ylabel('Cross-correlation');
ylim([-1 1])
title('Cross-correlation between average spike rate and eye trace during CS for suppressing neurons');
box('off');
set(gca, 'tickdir', 'out');
hold off

fname = sprintf('xcorr_spikes_eye_CS_suppression_%s.png', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
saveas(gca,file);
fname = sprintf('xcorr_spikes_eye_CS_suppression_%s.eps', group);
file = fullfile(p.figPath, "Cross-correlation",fname);
print(gcf, '-depsc', '-painters', file)


plotCorrelation(correlationData, "cs", group)
% plotCorrelation(correlationData, "us", group)

end

function trialFreq = calcAveSpkRateROI(data, roi)
p = IkUtils.getParams();
trialFreq = [];
for n2 = 1:length(data)
    rasterTrials = data(n2).RasterXY_cs_filtered(2,1:3:end);
    uniqueTrials = unique(sort(rasterTrials));
    trialFreqROI = [];
    for t = 1:220
        if ~ismember(t, uniqueTrials)
            trialFreqROI = [trialFreqROI;nan];
        else
            spikes = data(n2).RasterXY_cs_filtered(1,1:3:end);
            trialspikes = spikes(rasterTrials == t);
            range = roi.min:p.BinW:roi.max;
            trialspikecounts = histcounts(trialspikes, range);
            trialspikerate = sum(trialspikecounts)/(roi.max-roi.min);
            trialFreqROI = [trialFreqROI;trialspikerate];
        end
    end
    trialFreq = [trialFreq, {trialFreqROI}];
end

end

function [eye_trace_all, mean_spike_rate_roi] = calcStatsCC(correlationData, stim, range, mod)
if nargin < 4
    mod = "all";
end
p = IkUtils.getParams();
if stim == "cs"
    ROImin = p.tracesRanges.cs.min;
    ROImax = p.tracesRanges.cs.max;
elseif stim == "us"
    ROImin = p.tracesRanges.us.min;
    ROImax = p.tracesRanges.us.max;
end

mean_spike_rate_neuron = [];
eye_trace_all = [];

for cd = 1:length(correlationData.(stim).freq_trials)
    modMasksSes = correlationData.mod_masks{cd};
    cs_fac_mask = [];
    cs_sup_mask = [];
    us_fac_mask = [];
    us_sup_mask = [];
    for ms = 1:length(modMasksSes)
        cs_fac_mask = [cs_fac_mask, modMasksSes{ms}.cs_facilitation];
        cs_sup_mask = [cs_sup_mask, modMasksSes{ms}.cs_suppression];
        us_fac_mask = [us_fac_mask, modMasksSes{ms}.us_facilitation];
        us_sup_mask = [us_sup_mask, modMasksSes{ms}.us_suppression];
    end

    mean_eye_trace(cd,:) = correlationData(1).(stim).mean_eye_traces(cd,:);
    mean_eye_trace_roi(cd, :) = mean_eye_trace(cd, ROImin:ROImax);

    upsample_factor = 10;
    original_time = linspace(1, size(mean_eye_trace_roi(cd,:),2), size(mean_eye_trace_roi(cd,:),2));
    new_time = linspace(1, size(mean_eye_trace_roi(cd,:),2), size(mean_eye_trace_roi(cd,:),2)*upsample_factor);
    mean_eye_trace_upsampled(cd,:) = interp1(original_time, mean_eye_trace_roi(cd,:).', new_time, 'spline'); % Interpolate the data

    mean_spike_rate(cd) = [correlationData(1).mean_spike_rate(cd,:)];
    mean_spike_rate_neuron_ = mean_spike_rate{cd};

    bins = -p.t_raster_edges : p.spikeTimesBinW : p.t_raster_edges;
    
    rangeBase = find(bins < 0);
    test = correlationData.mean_spike_rate{cd};
    testBase = test(:,rangeBase);
    meanBase = mean(testBase,2);

    if mod == "fac" && stim == "cs"
        if ~isempty(mean_spike_rate_neuron_(logical(cs_fac_mask), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_(logical(cs_fac_mask), :) - meanBase(logical(cs_fac_mask))];
        end
        nModNeurons = sum(cs_fac_mask);
    elseif mod == "sup" && stim =="cs"
        if ~isempty(mean_spike_rate_neuron_(logical(cs_sup_mask), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_(logical(cs_sup_mask), :)- meanBase(logical(cs_sup_mask))];
        end
        nModNeurons = sum(cs_sup_mask);
    elseif mod == "fac" && stim == "us"
        if ~isempty(mean_spike_rate_neuron_(logical(us_fac_mask), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_(logical(us_fac_mask), :)- meanBase(logical(us_fac_mask))];
        end
        nModNeurons = sum(us_fac_mask);
    elseif mod == "sup" && stim =="us"
        if ~isempty(mean_spike_rate_neuron_(logical(us_sup_mask), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_(logical(us_sup_mask), :)- meanBase(logical(us_sup_mask))];
        end
        nModNeurons = sum(us_sup_mask);
    elseif mod == "nonmod" && stim =="cs"
        if ~isempty(mean_spike_rate_neuron_((~logical(cs_sup_mask) & ~logical(cs_fac_mask)), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_((~logical(cs_sup_mask) & ~logical(cs_fac_mask)), :)- meanBase((~logical(cs_sup_mask) & ~logical(cs_fac_mask)))];
        end
        nModNeurons = sum(~cs_fac_mask & ~cs_sup_mask);
    elseif mod == "nonmod" && stim =="us"
        if ~isempty(mean_spike_rate_neuron_((~logical(us_sup_mask) & ~logical(us_fac_mask)), :))
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_((~logical(us_sup_mask) & ~logical(us_fac_mask)), :)- meanBase((~logical(us_sup_mask) & ~logical(us_fac_mask)))];
        end
        nModNeurons = sum(~us_fac_mask & ~us_sup_mask);
    elseif mod == "all"
        if ~isempty(mean_spike_rate_neuron_)
        mean_spike_rate_neuron = [mean_spike_rate_neuron; mean_spike_rate_neuron_- meanBase];
        end
        nModNeurons = size(mean_spike_rate_neuron_,1);
    end

    for n3 = 1:nModNeurons %size(mean_spike_rate_neuron_,1)
        eye_trace_all = [eye_trace_all; mean_eye_trace_upsampled(cd,:)];
    end
end
try
mean_spike_rate_roi = mean_spike_rate_neuron(:,range);
catch
mean_spike_rate_roi = mean_spike_rate_neuron;
end
end

function plotCorrelation(correlationData, stim, group)
p = IkUtils.getParams();
% bins = -p.t_raster_edges : p.spikeTimesBinW : p.t_raster_edges;

% rangeBase = find(bins < 0);


if stim == "cs"
    ROImin = p.tracesRanges.cs.min;
    ROImax = p.tracesRanges.cs.max;
elseif stim == "us"
    ROImin = p.tracesRanges.us.min;
    ROImax = p.tracesRanges.us.max;
end
baseMin = p.tracesRanges.baseline.min;
baseMax = p.tracesRanges.baseline.max;
USmin = p.tracesRanges.us.min;
USmax = p.tracesRanges.us.max;

eyetraceRange = ROImin:ROImax;

p = IkUtils.getParams();

bins = -p.t_raster_edges : p.spikeTimesBinW : p.t_raster_edges;

rangeBase = find(bins < 0);
% mean_rates= correlationData.mean_spike_rate{1};
% mean_rates_base = mean_rates(:,rangeBase);
% meanBase = mean(mean_rates_base,2);

% Facilitation
fac_mask = [];
freq_fac_trials = [correlationData.(stim).freq_fac_trials];
for fft = 1: length(correlationData.(stim).freq_fac_trials)
    if ~isempty(freq_fac_trials{fft})
        fac_mask = [fac_mask; 1];
    else
        fac_mask = [fac_mask; 0];
    end
end
fac_mask = logical(fac_mask);
fac_idcs = find(fac_mask);

stimchar = char(stim);
stim1 = stimchar(1);

if ~isempty(fac_idcs)
    % Plot the scatter plot
    figure; hold on
    xlabel('Spike facilitation (Hz)');
    xlim([0 250])
    ylim([0 100])
    ylabel(sprintf('%sr amplitude percentage', stim1));
    title(sprintf('Rate and %sr amp %s facilitation', stim1, stim));
    title(sprintf('Plot of rate and %sr amplitude %s facilitation', stim1, stim));
box('off');
set(gca, 'tickdir', 'out');
    freq = [];
    CR = [];
    neuron_num = [];
    session_num = [];
    for fn = 1:length(fac_idcs)
        modMasksSes = correlationData.mod_masks{fac_idcs(fn)};
        cs_fac_mask = [];
        cs_sup_mask = [];
        us_fac_mask = [];
        us_sup_mask = [];
        for ms = 1:length(modMasksSes)
            cs_fac_mask = [cs_fac_mask, modMasksSes{ms}.cs_facilitation];
            cs_sup_mask = [cs_sup_mask, modMasksSes{ms}.cs_suppression];
            us_fac_mask = [us_fac_mask, modMasksSes{ms}.us_facilitation];
            us_sup_mask = [us_sup_mask, modMasksSes{ms}.us_suppression];
        end
        cs_fac_mask_idcs = find(cs_fac_mask);

        freq_neuron_ = freq_fac_trials{fac_idcs(fn)};
        eye_traces_neuron = correlationData.eye_traces{fac_idcs(fn)};
        CRamp_ROI = max(eye_traces_neuron(:, eyetraceRange)');
        baseline_min = nanmin(nanmean(eye_traces_neuron(:,baseMin:baseMax)));
        meanUR = nanmean(nanmax(eye_traces_neuron(:,USmin:USmax)'));
        fullBlinkRange = meanUR - baseline_min;
        CRamp_ROI = (CRamp_ROI-baseline_min) / fullBlinkRange * 100;

        CRamp_ROI(CRamp_ROI < p.thresCRperc) = nan;

        baserates = correlationData.mean_spike_rate{fac_idcs(fn)};

        for n = 1:size(freq_neuron_,2)
            meanBase = mean(baserates(n,rangeBase));

            freq_neuron = freq_neuron_{n};
            freq_neuron = freq_neuron';
            
            freq_neuron = freq_neuron - meanBase;

            mask = ~isnan(freq_neuron) & freq_neuron > 0 & ~isnan(CRamp_ROI);
            if sum(mask) > 1
            [R, P] = corrcoef(freq_neuron(mask), CRamp_ROI(mask));
            coefficients = polyfit(freq_neuron(mask), CRamp_ROI(mask), 1); % 1 indicates a linear fit
            y_fit = polyval(coefficients, freq_neuron(mask)); % Compute the corresponding y values
            alpha = 0.05; % Significance level
            if P(1,2) <= alpha
                line_color = 'r'; % Normal red for significant
                freq = [freq, {freq_neuron(mask)}];
                CR = [CR, {CRamp_ROI(mask)}];
                session_num = [session_num, fac_idcs(fn)];
                neuron_num = [neuron_num, cs_fac_mask_idcs(n)];
            else
                line_color = [1, 0.6, 0.6]; % Light red for non-significant
            end
                plot(freq_neuron(mask), y_fit, '-r', 'LineWidth', 2, 'Color', line_color);
            end
        end
    end
    hold off

    fname = sprintf('corr_linearReg_meanratetrials_maxamp_facilitating_%s_%s.png', stim, group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    saveas(gca,file);
    fname = sprintf('corr_linearReg_meanratetrials_maxamp_facilitating_%s_%s.eps', stim, group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    print(gcf, '-depsc', '-painters', file)

for f = 1 : numel(freq)
    figure();
    hold on
    freqLoop = freq{f};
    CRLoop = CR{f};
    [R, P] = corrcoef(freqLoop, CRLoop);
box('off');
set(gca, 'tickdir', 'out');
    scatter(freqLoop, CRLoop, 'filled'); % Use 'filled' to fill the markers
    coefficients = polyfit(freqLoop, CRLoop, 1); % 1 indicates a linear fit
    % Evaluate the fit
    y_fit = polyval(coefficients, freqLoop); % Compute the corresponding y values
    plot(freqLoop, y_fit, '-r', 'LineWidth', 2);

    % Display the correlation coefficient on the plot
    
    xlabel('frequency');
    ylabel(sprintf('%sr amplitude', stim1));
    title(sprintf('Scatter Plot of rate and %sr amplitude %s facilitation', stim1, stim));
    try
        var = ['Correlation: ', num2str(R(1,2)), '   P: ', num2str(P(1,2))];
        text(min(freqLoop), max(CRLoop), var, 'FontSize', 12, 'Color', 'red');
    catch
    end

    fname = sprintf('correlation_single_neuron_facilitating_%s_S%d_N%d.png', group, session_num(f), neuron_num(f));
    file = fullfile(p.figPath, "Cross-correlation",fname);
    saveas(gca,file);
    fname = sprintf('correlation_single_neuron_facilitating_%s_S%d_N%d.eps', group, session_num(f), neuron_num(f));
    file = fullfile(p.figPath, "Cross-correlation",fname);
    print(gcf, '-depsc', '-painters', file)

end

end


%% suppressing

sup_mask = [];
freq_sup_trials = [correlationData.(stim).freq_sup_trials];
for fft = 1: length(correlationData.(stim).freq_sup_trials)
    if ~isempty(freq_sup_trials{fft})
        sup_mask = [sup_mask; 1];
    else
        sup_mask = [sup_mask; 0];
    end
end
sup_mask = logical(sup_mask);
sup_idcs = find(sup_mask);


freq = [];
CR = [];
neuron_num = [];
session_num = [];
if ~isempty(sup_idcs)
    % Plot the scatter plot
    figure; hold on
    xlabel('Spike suppression (Hz)');
    xlim([-150 0])
    ylim([0 100])
    ylabel(sprintf('%sr amplitude percentage', stim1));
    title(sprintf('Rate and %sr amp %s suppression', stim1, stim));
box('off');
set(gca, 'tickdir', 'out');
    for fn = 1:length(sup_idcs)
        modMasksSes = correlationData.mod_masks{sup_idcs(fn)};
        cs_fac_mask = [];
        cs_sup_mask = [];
        us_fac_mask = [];
        us_sup_mask = [];
        for ms = 1:length(modMasksSes)
            cs_fac_mask = [cs_fac_mask, modMasksSes{ms}.cs_facilitation];
            cs_sup_mask = [cs_sup_mask, modMasksSes{ms}.cs_suppression];
            us_fac_mask = [us_fac_mask, modMasksSes{ms}.us_facilitation];
            us_sup_mask = [us_sup_mask, modMasksSes{ms}.us_suppression];
        end
        cs_sup_mask_idcs = find(cs_sup_mask);

        freq_neuron_ = freq_sup_trials{sup_idcs(fn)};
        eye_traces_neuron = correlationData.eye_traces{sup_idcs(fn)};
        CRamp_ROI = max(eye_traces_neuron(:, eyetraceRange)');
        baseline_min = nanmin(nanmean(eye_traces_neuron(:,baseMin:baseMax)));
        meanUR = nanmean(nanmax(eye_traces_neuron(:,USmin:USmax)'));
        fullBlinkRange = meanUR - baseline_min;
        CRamp_ROI = (CRamp_ROI-baseline_min) / fullBlinkRange * 100;

        CRamp_ROI(CRamp_ROI < p.thresCRperc) = nan;
        baserates = correlationData.mean_spike_rate{sup_idcs(fn)};

        for n = 1:size(freq_neuron_,2)
            meanBase = mean(baserates(n,rangeBase));
            freq_neuron = freq_neuron_{n};
            freq_neuron = freq_neuron';
            freq_neuron = freq_neuron - meanBase;

            mask = ~isnan(freq_neuron) & freq_neuron < 0 & ~isnan(CRamp_ROI);
            if sum(mask) > 1
            [R, P] = corrcoef(freq_neuron(mask), CRamp_ROI(mask));
            coefficients = polyfit(freq_neuron(mask), CRamp_ROI(mask), 1); % 1 indicates a linear fit
            y_fit = polyval(coefficients, freq_neuron(mask)); % Compute the corresponding y values
            alpha = 0.05; % Significance level
            if P(1,2) <= alpha
                line_color = 'r'; % Normal red for significant
                freq = [freq, {freq_neuron(mask)}];
                CR = [CR, {CRamp_ROI(mask)}];
                session_num = [session_num, sup_idcs(fn)];
                neuron_num = [neuron_num, cs_sup_mask_idcs(n)];        
            else
                line_color = [1, 0.6, 0.6]; % Light red for non-significant
            end
            plot(freq_neuron(mask), y_fit, '-r', 'LineWidth', 2, 'Color', line_color);
            end
        end
    end
    hold off

    fname = sprintf('corr_linearReg_meanratetrials_maxamp_suppressing_%s_%s.png',stim, group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    saveas(gca,file);
    fname = sprintf('corr_linearReg_meanratetrials_maxamp_suppressing_%s_%s.eps',stim, group);
    file = fullfile(p.figPath, "Cross-correlation",fname);
    print(gcf, '-depsc', '-painters', file)

end

for f = 1 : numel(freq)
    figure();
    hold on
    freqLoop = freq{f};
    CRLoop = CR{f};
    [R, P] = corrcoef(freqLoop, CRLoop);
box('off');
set(gca, 'tickdir', 'out');
    scatter(freqLoop, CRLoop, 'filled'); % Use 'filled' to fill the markers
    coefficients = polyfit(freqLoop, CRLoop, 1); % 1 indicates a linear fit
    % Evaluate the fit
    y_fit = polyval(coefficients, freqLoop); % Compute the corresponding y values
    plot(freqLoop, y_fit, '-r', 'LineWidth', 2);

    % Display the correlation coefficient on the plot
    
    xlabel('frequency');
    ylabel(sprintf('%sr amplitude', stim1));
    title(sprintf('Scatter Plot of rate and %sr amplitude %s suppression', stim1, stim));
    try
        var = ['Correlation: ', num2str(R(1,2)), '   P: ', num2str(P(1,2))];
        text(min(freqLoop), max(CRLoop), var, 'FontSize', 12, 'Color', 'red');
    catch
    end

    fname = sprintf('correlation_single_neuron_suppressing_%s_S%d_N%d.png', group, session_num(f), neuron_num(f));
    file = fullfile(p.figPath, "Cross-correlation",fname);
    saveas(gca,file);
    fname = sprintf('correlation_single_neuron_suppressing_%s_S%d_N%d.eps', group, session_num(f), neuron_num(f));
    file = fullfile(p.figPath, "Cross-correlation",fname);
    print(gcf, '-depsc', '-painters', file)

end


end
