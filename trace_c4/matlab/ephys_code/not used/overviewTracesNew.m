clear all
close all

%% preallocation
p = IkUtils.getParams();
base_dir = p.pathBehaviorDataHome;
lastDay = 2;
nTimeSteps = 200;
nTrials = 240;
nCSplusRegTrials = 220;

defaultMiceList = defaultMice();
codes = [defaultMiceList([defaultMice().type] == "Shank2MUT")];
codes = [codes.code];
MUTmice = codes;
MUTmice = MUTmice(1:end-1);
codes = [defaultMiceList([defaultMice().type] == "Shank2WT")];
codes = [codes.code];
WTmice = codes;
WTmice = WTmice(1:end-1);

type = ones(1,nTrials);

WT = struct('traces', cell(1, numel(WTmice)), ...
            'cramp', cell(1, numel(WTmice)), ...
            'cramp5', cell(1, numel(WTmice)), ...
            'mean', cell(1, numel(WTmice)), ...
            'onsets', cell(1, numel(WTmice)));

MUT = struct('traces', cell(1, numel(MUTmice)), ...
             'cramp', cell(1, numel(MUTmice)), ...
             'cramp5', cell(1, numel(MUTmice)), ...
             'mean', cell(1, numel(MUTmice)), ...
             'onsets', cell(1, numel(MUTmice)));

%% calculations - WT
for i = 1:numel(WTmice)
    WT(i).traces = nan(lastDay, nTrials * nTimeSteps);
    WT(i).cramp = nan(lastDay, nTrials, nCSplusRegTrials);
    WT(i).cramp5 = nan(lastDay, nTrials, nCSplusRegTrials);
    WT(i).mean = nan(lastDay, nTimeSteps);
    WT(i).onsets = nan(lastDay, nTrials);
    
    mouse = WTmice(i);
    dir_content = dir(fullfile(base_dir, mouse));
    path_list = IkUtils.getPathFromDir(dir_content);
    
    for f = 1:numel(path_list)
        folder = path_list(f);
        load(fullfile(folder, 'trialdata.mat'));
        WT(i).traces(f, :) = reshape(behavior_trial_data.tracesnorm(1:nTrials, :), 1, []);
        WT(i).cramp(f, :, :) = permute(behavior_trial_data.CRamp(1:nCSplusRegTrials, :), [3, 1, 2]);
        WT(i).cramp5(f, :, :) = permute(behavior_trial_data.CRamp5, [3, 1, 2]);
        
        WT(i).onsets(f, :) = behavior_trial_data.tm(1, find(any(behavior_trial_data.tracesnorm > behavior_trial_data.tracesnorm(:, 1:40) + 5 * std(behavior_trial_data.tracesnorm(:, 1:40)), 2), 1));
    end
end

%% calculations - MUT
for i = 1:numel(MUTmice)
    MUT(i).traces = nan(lastDay, nTrials * nTimeSteps);
    MUT(i).cramp = nan(lastDay, nTrials, nCSplusRegTrials);
    MUT(i).cramp5 = nan(lastDay, nTrials, nCSplusRegTrials);
    MUT(i).mean = nan(lastDay, nTimeSteps);
    MUT(i).onsets = nan(lastDay, nTrials);
    
    mouse = MUTmice(i);
    dir_content = dir(fullfile(base_dir, mouse));
    path_list = IkUtils.getPathFromDir(dir_content);
    
    for f = 1:numel(path_list)
        folder = path_list(f);
        load(fullfile(folder, 'trialdata.mat'));
        MUT(i).traces(f, :) = reshape(behavior_trial_data.tracesnorm(1:nTrials, :), 1, []);
        MUT(i).cramp(f, :, :) = permute(behavior_trial_data.CRamp(1:nCSplusRegTrials, :), [3, 1, 2]);
        MUT(i).cramp5(f, :, :) = permute(behavior_trial_data.CRamp5, [3, 1, 2]);
        
        MUT(i).onsets(f, :) = behavior_trial_data.tm(1, find(any(behavior_trial_data.tracesnorm > behavior_trial_data.tracesnorm(:, 1:40) + 10 * std(behavior_trial_data.tracesnorm(:, 1:40)), 2), 1));
    end
end

%% Final figures
% You can keep your existing plotting code here.

print(eyeblinkTracesOverview, fullfile(base_dir, 'eyeblinkTracesOverview.png'), '-dpng')
print(CROverview, fullfile(base_dir, 'CROverview.png'), '-dpng')
