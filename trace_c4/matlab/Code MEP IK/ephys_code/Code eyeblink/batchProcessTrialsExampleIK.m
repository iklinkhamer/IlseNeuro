%% IK 2-11-23
function batchProcessTrialsExampleIK(mInput)
close all
app = getUImain(scriptname = mfilename);

P = IkUtils.getParams();

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    disp("Calculating...")
    appData = evalin('base', 'appData');

    startMouseCode = appData.startMouseCode;
    onlySeeFigs = appData.seeFigs;
    overrideFigs = appData.overrideFigs;
    debug = appData.debugTraces;
    plotfiguresOld = appData.plotFigsOld;
    onlySaveFigs = appData.onlySaveFigs;

    mice = defaultMice();
    mouseNames = [mice.name];
    mouseCodeList = [mice.code];
    if nargin < 1
        mouseCodes = appData.mousecodes;
    else
        mouseCodes = IkUtils.promptMouseNames(mInput);
    end
    startMouseIdx = find(mouseCodes == startMouseCode);

    if length(mouseCodes) == 1
        startMouseIdx = 1;
    end

    for m = startMouseIdx:length(mouseCodes)
        close all
        mcode = mouseCodes(m);
        path_list = getBehaviorPathList(mcode);
        figpath = fullfile(P.figPath, "Training per mouse", mcode);

        for f = 1 : length(path_list)
            folder = path_list(f);
            folder_contents = dir(fullfile(folder, 'trialdata*.mat'));
            [~, folderName, ~] = fileparts(folder);
            if contains(folderName, P.s)
                s = find(P.s == folderName);
            else
                s = folderName;
            end

            if (onlySaveFigs || onlySeeFigs) || debug == 1
                go = 1;
            else
                continue
            end
            fprintf("Mouse: %s, %s Session: %s\n", mcode, mouseNames(mouseCodeList == mcode), string(s));
            stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');
            if isempty(folder_contents) || debug == 1 && ~(isempty(dir(fullfile(folder, "*encoder.mat"))))
                behavior_trial_data = processTrials(folder, 'recalibrate'); % Recalibrate eyelid
            elseif go
                try
                    load(fullfile(folder, stamptrial));
                catch
                    try
                        load(fullfile(folder, 'trialdata.mat'))
                    catch
                        try
                        catch
                            load(fullfile(folder, folder_contents(1).name));
                            continue
                        end
                    end
                end
            else
                continue
            end


            if ~isempty(behavior_trial_data)

                results_folder = getResultsFolder(folder, mcode);
                try
                    [eyelidpos, eye_closing_mask, eye_opening_mask] = upsampleEyelidpos(behavior_trial_data);
                catch
                    disp("Probably something wrong with eye traces file.")
                    continue
                end
                stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyelidpos', file_extension = '.mat');
                save(fullfile(folder, stampeye), 'eyelidpos');
                stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyeClosingMask', file_extension = '.mat');
                save(fullfile(folder, stampeye), 'eye_closing_mask');
                stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyeOpeningMask', file_extension = '.mat');
                save(fullfile(folder, stampeye), 'eye_opening_mask');
                if ~exist(fullfile(results_folder, stampeye), 'file')
                    save(fullfile(results_folder, stampeye), 'eyelidpos');
                end

                stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');
                if ~exist(fullfile(results_folder, stamptrial), 'file')
                    save(fullfile(results_folder, stamptrial), 'behavior_trial_data');
                end

                %             plotFigureEyeTraces(behavior_trial_data.tm(1,:), smoothedTraces);
                %                                     plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces)

                if ~exist(figpath, 'dir')
                    mkdir(figpath)
                end
                if overrideFigs || ~exist(fullfile(figpath, sprintf("EyelidClosure_%s_%s.jpg", mcode, folderName)), 'file')
                    traces = behavior_trial_data.eyelidpos;
                    smoothedTraces = smootheTraces(traces);
                    normalized_traces = normalizeEyelidTraces(smoothedTraces, behavior_trial_data);
                    fig = plotFigureEyeTraces(behavior_trial_data.tm(1,:), normalized_traces);
                    saveas(fig, fullfile(figpath, sprintf("EyelidClosure_%s_%s.jpg", mcode, folderName)));
                end
                if overrideFigs || ~exist(fullfile(figpath, sprintf("EyelidClosure_%s_%s.eps", mcode, folderName)), 'file')
                    traces = behavior_trial_data.eyelidpos;
                    smoothedTraces = smootheTraces(traces);
                    normalized_traces = normalizeEyelidTraces(smoothedTraces, behavior_trial_data);
                    fig = plotFigureEyeTraces(behavior_trial_data.tm(1,:), normalized_traces);
                    print(fig, '-depsc', '-painters', fullfile(figpath, sprintf("EyelidClosure_%s_%s.eps", mcode, folderName)))
                end
                if onlySeeFigs && ~onlySaveFigs
                    traces = behavior_trial_data.eyelidpos;
                    smoothedTraces = smootheTraces(traces);
                    normalized_traces = normalizeEyelidTraces(smoothedTraces, behavior_trial_data);
                    fig = plotFigureEyeTraces(behavior_trial_data.tm(1,:), normalized_traces);
                end


                if plotfiguresOld
                    plotFiguresOld(behavior_trial_data)
                end
            end
        end
    end
    disp("Loop end. Waiting for new UI input...")
end

end

function results_folder = getResultsFolder(folder, mcode)
[~, folderName, ~] = fileparts(folder);
% if ~regexp(folder, "TrainingOnly")
    resultsFolderPattern = "behaviorDataResults";
% else
%     resultsFolderPattern = "behaviorDataResultsTrainingOnlyMice";
% end
behaviorDataResultsFolder = fullfile("/home/i.klinkhamer/Documents/Data/", resultsFolderPattern, mcode, folderName);
if ~exist(behaviorDataResultsFolder, 'dir')
    mkdir(behaviorDataResultsFolder);
end

first_part_folder = fileparts(fileparts(fileparts(folder)));
middle_folder = resultsFolderPattern;
end_folder = erase(folder, fileparts(fileparts(folder)));
results_folder = fullfile(first_part_folder,middle_folder,end_folder);
end

function [eyelidpos, eye_closing_mask, eye_opening_mask] = upsampleEyelidpos(behavior_trial_data)
eyelidpos = behavior_trial_data.eyelidpos;
% Define the factor by which you want to upsample
upsample_factor = 10;
% Create the original time points
original_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2));
% Create the new time points after upsampling
new_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2)*upsample_factor);
% Interpolate the data
eyelidpos_upsampled = interp1(original_time, eyelidpos.', new_time, 'spline');
smoothedTraces = smootheTraces(eyelidpos_upsampled);
smoothedTraces = smoothedTraces';
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max*10;
USmin = baseMax + 1;
USmax = P.tracesRanges.us.max*10;

baseline_min = min(nanmean(smoothedTraces(:, baseMin:baseMax), 2));      % Calculate baselines as the mean of traces within the specified range
traces_2baseline = smoothedTraces - ...
    (nanmean(smoothedTraces(:,baseMin:baseMax), 2) - ...
    baseline_min);
min_val = min(nanmean(traces_2baseline(:,baseMin:baseMax)'));
max_val = max(max(traces_2baseline(:,USmin:USmax)'));
normalized_traces = (traces_2baseline - min_val) / (max_val - min_val);
eyelidpos = normalized_traces;

dy_dx = gradient(eyelidpos);
window_size = 10;
avg_gradient = movmean(dy_dx, window_size);

threshold = 0.001;
% Determine if the gradient is positive or negative based on the average
eye_closing_mask = avg_gradient > threshold;
eye_opening_mask = avg_gradient < -threshold;

eyelidpos = eyelidpos';
eye_closing_mask = eye_closing_mask';
eye_opening_mask = eye_opening_mask';
end

function fig = plotFigureEyeTraces(tm, traces)
fig = figure; hold on
P = IkUtils.getParams();
x1cs = P.eventTimes(1);
x2cs = P.eventTimes(2);
y1 = -0.1;
y2 = 1;
num_lines = length(traces);
cmap = copper(num_lines);
mintrans = 0.4;
transparancy = flip(mintrans:(1-mintrans)/(num_lines+1):1);
patch ...
    ( ...
    [x1cs x1cs x2cs x2cs] ...
    , [y1 y2 y2 y1] ...
    , 'blue' ...
    , EdgeColor = 'none' ...
    , FaceAlpha = 0.2 ...
    )
x1us = P.eventTimes(2);
x2us = x1us + P.USdur;
patch ...
    ( ...
    [x1us x1us x2us x2us] ...
    , [y1 y2 y2 y1] ...
    , 'green' ...
    , EdgeColor = 'none' ...
    , FaceAlpha = 0.2 ...
    )
for t = 1 : size(traces,1)
    plot(tm, traces(t,:), 'Color', [cmap(t,:), transparancy(t)])
end
meanTrace = nanmean(traces(:,:));
plot(tm, meanTrace, 'Color', 'k', 'LineWidth', 2)
xlabel("Time (ms)")
ylabel("Eyelid closure ratio")
ylim([y1 y2])
yticks(0:0.2:y2)
% yticklabels()
set(gca, "TickDir", 'out')
box('off')
cbar = colorbar; % Create colorbar
cbar.Ticks = linspace(0, 40, num_lines+1); % Set colorbar ticks
cbar.TickLabels = 0:40:num_lines; % Set colorbar tick labels
cbar.Colormap = cmap; % Set colorbar colormap to match line colors
ylabel(cbar, 'Trial number');
title("Eyelid closure")

end

function plotFiguresOld(behavior_trial_data)
baseline = mean(mean(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,1:40)));
traces_new = behavior_trial_data.eyelidpos - (mean(behavior_trial_data.eyelidpos(:,1:40), 2) - baseline);
baseline_new = min(mean(traces_new(behavior_trial_data.c_csdur==0,1:40),2));
stdev = std(behavior_trial_data.eyelidpos(:,1:40), 0, 2);
mask_stdev = stdev < 3*mean(stdev);
mask_lowest = min(traces_new(:,41:end), [], 2) > (baseline_new - 0.3);
mask_highest = max(traces_new(:,141:end), [],2) < max(max(traces_new(:,90:140)'));
traces_new = traces_new(mask_stdev & mask_lowest & mask_highest, :);
min_val = min(traces_new(:)');
max_val = max(traces_new(:)');
traces_new = (traces_new - min_val) / (max_val - min_val);
plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces_new)
end
