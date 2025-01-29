%% IK 8-4-24
load('/home/i.klinkhamer/Documents/Data/spikeSortingUnits/Shank2MUT/S1/trialspikes.mat')
p = IkUtils.getParams();
n_neurons = length(trialspikes);
n_trials = length(trialspikes(1).spiketimes_cs_neuron);
bins = p.psthRanges.cs_full.min:p.spikeTimesBinW:p.psthRanges.cs_full.max;
spikecounts = zeros(length(bins)-1, n_neurons*n_trials);
i = 1;
for n = 1 : n_neurons
    for t = 1 : n_trials
        spikes = trialspikes(n).spiketimes_cs_neuron(t).trial;
        spikecounts(:,i) = histcounts(spikes,bins).';
        i = i + 1;
    end
end
save(fullfile("/home/i.klinkhamer/Documents/Data/spikeSortingUnits/Shank2MUT/S1/",'spikecountsallneuronsalltrials.mat'),'spikecounts')
%%
p = IkUtils.getParams();

mouseList = p.mouseList();

for i = 1:length(mouseList)
    status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s", mouseList(i)));
    status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S1", mouseList(i)));
    status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S2", mouseList(i)));
    if i > 52
        status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S3", mouseList(i)));
    end
end
%
for i = 1:length(mouseList)
    if isfolder(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1", mouseList(i)))
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/RAWdata.bin", mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/settings.xml",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/stimtimes.mat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/timestamps.mat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/continuous.dat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
    end
    if isfolder(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2", mouseList(i)))
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/RAWdata.bin", mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/settings.xml",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/stimtimes.mat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/timestamps.mat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
        status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/continuous.dat",...
            mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
    end
end
%%
sessions = [1, 2];
for i = 1:length(mouseList)
    if isfolder(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s", mouseList(i)))
        %         status = mkdir(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/S1", mouseList(i)));
        a = dir(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/", mouseList(i)));
        for j = 1:length(a)
            status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_1/%s", mouseList(i), a(j).name),...
                sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
        end
    end
    if isfolder(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s", mouseList(i)))
        %         status = mkdir(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/S2", mouseList(i)));
        a = dir(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/", mouseList(i)));
        for j = 1:length(a)
            status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/Session_2/%s", mouseList(i), a(j).name),...
                sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
        end
    end
end
%%
for i = 1:length(mouseList)
    status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/01", mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S1", mouseList(i)));
    status = movefile(sprintf("/home/mick/Desktop/Ilse/spikeSortingUnits/%s/02", mouseList(i)), sprintf("/home/mick/Desktop/Ilse/convertedData/%s/S2", mouseList(i)));
end


figure;
ButtonH = uicontrol('Style', 'ToggleButton', 'String', 'Stop', ...
    'Value', 0);
for i = 1:1e6
    drawnow;
    if get(ButtonH, 'Value') == 1
        keyboard;
    end
end

keyboard_prompt = input("Give control to keyboard? ");
if keyboard_prompt == 1
    keyboard
end
%%
spikes = [];

spikes = [spikes; unit.neuron(n).trial_spikes(trial).trial];

spikes = sort(spikes);
idx_cs_new = find(spikes > -0.15 & spikes < -0.05);

RasterXY_cs_CSUS = unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CSUS_trials));
RasterXY_cs_CS_only = unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CS_only_trials));
idx_cs_CSUS = find(RasterXY_cs_CSUS > -0.15 & RasterXY_cs_CSUS < -0.05);
idx_cs_CS_only = find(RasterXY_cs_CS_only > -0.15 & RasterXY_cs_CS_only < -0.05);

fr_CSUS = filter(b,a,unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CSUS_trials)));
fr_CS = filter(b,a,unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CS_only_trials)));
unit.neuron(n).lines_cs_reset = zeros(1,(length(CSUS_trials) + length(CS_only_trials)));
unit.neuron(n).lines_cs_reset(CSUS_trials) = fr_CSUS - mean(fr_CSUS(idx_cs_CSUS)) + base_fr  ;
unit.neuron(n).lines_cs_reset(CS_only_trials) = fr_CS - mean(fr_CS(idx_cs_CS_only)) + base_fr  ;


plot(unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CSUS_trials)),...
    unit.neuron(n).lines_cs_reset(CSUS_trials),...
    unit.neuron(n).RasterXY_cs(1,ismember(unit.neuron(n).RasterXY_cs(2,:), CS_only_trials)),...
    unit.neuron(n).lines_cs_reset(CS_only_trials));


%% I.K. 8/9/23   PSTH for combined US/CS
function makeRasterPlots(kwargs)
arguments
    kwargs.loadDataFromDATA = true;
    kwargs.convertAllSessions = true;
    kwargs.convertAllMice = true;
    kwargs.startMouse = 51;
end
close all
addpath('/media/mick/DATA/Ilse/convertedData/')

P = getParams();

if kwargs.convertAllMice == 0
    prompt = "Enter mouse index: ";
    mousenames = string(input(prompt, 's'));
else
    mousenames = P.mouseList;
    mousenames = mousenames(kwargs.startMouse:end);
end

if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    sessions = 1:P.n_sessions;
end

for mname = mousenames
    for s = 1:length(sessions)
        fprintf("Mouse: %s Session: %d \n", mname, s);
        [unit, S] = getData(kwargs, mname, s);
        if isempty(unit)
            continue
        end

        stimtimes = unit.stimtimes_sorted; % IK note: uncomment this if
        %         you want to see the CS only trials at the bottom.
        %         stimtimes = unit.stimtimes;
        n_trials = length(stimtimes);
        n_neurons = length(unit.neuron);

        % extract spiketimes and stimtimes
        for n = 1 : n_neurons
            unit.neuron(n).trial_numbers = zeros(length(unit.neuron(n).spiketimes),1);
            for trial = 1 : n_trials
                unit.neuron(n).trial_numbers(unit.neuron(n).spiketimes > (stimtimes(trial) - P.t_pre_trial) &...
                    unit.neuron(n).spiketimes < (stimtimes(trial) + P.t_post_trial)) =...
                    trial;
            end
            mask_trials = unit.neuron(n).trial_numbers > 0;
            unit.neuron(n).RasterXY_cs = zeros(2, sum(mask_trials)*3);
            unit.neuron(n).RasterXY_cs(2, 1:3:end) = unit.neuron(n).trial_numbers(mask_trials);
            unit.neuron(n).RasterXY_cs(2, 2:3:end) = unit.neuron(n).trial_numbers(mask_trials) + P.raster_spike_height;
            unit.neuron(n).RasterXY_cs(2, 3:3:end) = nan;
            unit.neuron(n).RasterXY_cs(1, 3:3:end) = nan;
        end

        for n = 1 : n_neurons
            for trial = 1 : n_trials
                unit.neuron(n).trial_spikes(trial).trial = unit.neuron(n).spiketimes(unit.neuron(n).trial_numbers == trial) - stimtimes((trial));
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial) = unit.neuron(n).trial_spikes(trial).trial;
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial + P.raster_spike_height) = unit.neuron(n).trial_spikes(trial).trial;
            end
        end

        % psth figure

        for n = 1 : n_neurons
            psth = figure(1);
            set(psth, 'Color', 'white','Position', [0 2000 1800 1000]);
            subplot(1,3,1)
            hold on
            [Spike_counts,bins] = histcounts(unit.neuron(n).RasterXY_cs(1, 1:3:end),100);

            h = histogram('BinCounts', Spike_counts, 'BinEdges', bins, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0.9, 0.9, 0.9], 'FaceAlpha', 0.8);

            % Add presentation markers for CS and US
            a1 = area([0 0.27], [max(Spike_counts) max(Spike_counts)]);
            a1.FaceColor = [0,0,1];
            a1.FaceAlpha = 0.3;

            a2 = area([0.25 0.27], [max(Spike_counts) max(Spike_counts)]);
            a2.FaceColor = [0,1,0];
            a2.FaceAlpha = 0.3;

            % Plot lines to mark the presentation of CS and US
            hold on
            plot([0 0.27], [max(Spike_counts)*1.1 max(Spike_counts)*1.1], 'LineWidth', 2, 'Color', [0,0,1])
            plot([0.25 0.27], [max(Spike_counts)*1.1 max(Spike_counts)*1.1], 'LineWidth', 2, 'Color', [0,1,0])

            ymax = max(Spike_counts)*1.2;

            % Set axis limits and labels
            hold off
            xlim([-0.5 2])
            ylim([0 ymax])
            xlabel('Time (s)', 'FontSize', 14)
            ylabel('Spike Count', 'FontSize', 14)



            % Add a title to the plot
            title(sprintf('PSTH of neuron %d', n), 'FontSize', 16)
            if n ==1
                legend( 'CS presented', 'US presented', 'Location', 'southeast')
            end


            % CS- Aligned Raster
            %         raster = figure('Color', 'white','Position', [10 0 1500 700]);
            subplot(1,3,2)
            hold on

            title(sprintf('CS Aligned: Sess %d  N  %d', s, unit.irc_trueID(n)))
            plot(unit.neuron(n).RasterXY_cs(1,:),unit.neuron(n).RasterXY_cs(2,:),'Color',[0.2 0.2 0.2])
            mn = 1; mx = P.n_trials;
            plot([0 0], [mn mx],'k--')

            plot([P.t_US P.t_US], [mn mx],'k--')

            xlim([-0.5 2]); ylim([1 P.n_trials]);
            hold off

            % CS-aligned lines
            subplot(1,3,3), hold on
            [a, b] = GetButter();
            RasterXY_spikes = unit.neuron(n).RasterXY_cs(1,1:3:end);
            RasterXY_trials = unit.neuron(n).RasterXY_cs(2,1:3:end);
            idx_cs = find(unit.neuron(n).RasterXY_cs(1,1:3:end) > -0.15 & unit.neuron(n).RasterXY_cs(1,1:3:end) < -0.05);
            base_fr = mean(mean(RasterXY_spikes(1,idx_cs)));
            for c = 1 : P.n_conditions
                trials = find(stimtimes(:,2) == P.conditions(c));
                for trial = 1 : length(trials)
                    [Spike_counts,bins] = histcounts(RasterXY_spikes(ismember(RasterXY_trials, trials)),100);
                end
                mask_cs_bins = bins > -0.15 & bins < -0.05;
                fr = filter(b,a,Spike_counts);
                unit.neuron(n).lines_cs_reset(c).fr = fr - mean(fr(mask_cs_bins)) + base_fr  ;

                Spike_counts = Spike_counts./length(trials);
                % create filter
                sigma = 2; % pick sigma value for the gaussian
                gaussFilter = gausswin(6*sigma + 1)';
                gaussFilter = gaussFilter / sum(gaussFilter); % normalize
                filteredY = conv(Spike_counts, gaussFilter, 'same');


                % plot(bins(1:end-1), Spike_counts)
                plot(bins(1:end-1), filteredY)
            end
            ymax = max(Spike_counts)*1.2;
            mn = 0; mx = ymax;
            plot([0 0], [mn mx],'k--')

            plot([P.t_US P.t_US], [mn mx],'k--')

            xlim([-0.5 2]); ylim([0 ymax]);
            hold off

            pause

            figure(1),subplot(1,2,1), cla; figure(1),subplot(1,2,2), cla;
        end
    end
end


%% saving

% save(strcat(config_path,'spikesCSUS.mat'),'Spikes')
%print(psth, fullfile(config_path, 'psth.png'), '-dpng')


end

function [unit, S] = getData(kwargs, mname, s)
P = getParams();
directory = P.dirHome;
path_end = mname + "/" + P.s(s) + "/";
config_path = directory + path_end;
if ~isfolder(config_path) && s == 2
    unit = [];    S = [];
    disp("Only one session found")
    return
end
filename = getFilename(config_path);
if isempty(filename)
    unit = [];    S = [];
    disp("csv file not found")
    return
end
if kwargs.loadDataFromDATA == 1
    dataPath = P.dataPathDATA + path_end;
else
    dataPath = P.dataPathHome + path_end;
end
if ~isfile(strcat(dataPath,'stimtimes.mat'))
    unit = []; S = [];
    disp("No stimtimes found, this is a Rhymdata mouse.")
    return
end
cell_spk = import_jrc_csv(strcat(filename));
n_neurons = length(cell_spk) - 1;
if n_neurons == 0
    unit = []; S = [];
    disp("There are no neurons in this session")
    return
end

S = load(strcat(erase(filename, ".csv"),"_jrc.mat"));

for i = 1 : n_neurons
    unit.neuron(i).spiketimes = cell_spk(i+1).t;
    unit.neuron(i).irc_trueID = cell_spk(i+1).nr;
    unit.neuron(i).clusters = cell_spk(i+1).ch;
    unit.irc_trueID(i) = cell_spk(i+1).nr;
end

unit.spiketimes = cell_spk(1).t;


stimtimes = load(strcat(dataPath,'stimtimes.mat'));
stimtimes = stimtimes.times(stimtimes.times <= 5000); % IK note: 3500
% stimtimes = [stimtimes zeros(length(stimtimes),1)];
stimtimes = stimtimes(1:2:end); % Take only the odd stimtimes, because those are the times that the light is turning ON.
try
    stimtimes = stimtimes(1:480, :);
catch
    try
        stimtimes = stimtimes(1:440, :);
    catch
    end
end
%stimtimes =  stimtimes(1:11:end,:);
% stimtimes(1:11:end,:) = []; % IK note: uncomment this if you only want
% the control CS-only trials
% stimtimes(1:2:end,2) = 1;
%         stimtimes = stimtimes(1:2:end,:);
% stimtimes(1:11:end,3) = 1; % IK: CS only trials.
% stimtimes(1:12:end,3) = 1; % IK: CS only trials.
stimtimes(:,2) = 1;
stimtimes(1:5:end,2) = 0; % IK: CS only trials.
% stimtimes_CS_only_trials_first = sortrows(stimtimes, 3, 'descend');
stimtimes_CS_only_trials_first = sortrows(stimtimes, 2, 'ascend');


unit.stimtimes = stimtimes;
unit.stimtimes_sorted = stimtimes_CS_only_trials_first;
end

function filename = getFilename(csv_path)
if isfile(csv_path + "Config_h32_oe_bin.csv")
    filename = csv_path + "Config_h32_oe_bin.csv";
elseif isfile(csv_path + "Config_h32_oe_dat.csv")
    filename = csv_path + "Config_h32_oe_dat.csv";
elseif isfile(csv_path + "Config_h32_oe.csv")
    filename = csv_path + "Config_h32_oe.csv";
else
    filename = [];
end
end


% convertData

dataPathDATA = sprintf('%s/%s/%s/', P.pathEphysDataDATA, mname, P.s(s));

path = sprintf('/mnt/Data/Ilse/Data/%s/', mname);
dataPathHome = sprintf('%s/%s/%s/', P.pathEphysDataHome, mname, P.s(s));

if isfile(dataPathHome + "all_channels.events") || isfile(dataPathHome + "all_channels_2.events")
    path = dataPathHome;
else
    path = dataPathDATA;
end

% initSpikeSorting

dataPathDATA = P.pathEphysDataDATA + mname + "/" + session_n + "/";
dataPathHome = P.pathEphysDataHome + mname + "/" + session_n + "/";
if isfile(dataPathHome + "RAWdata.bin") || isfile(dataPathHome + "continuous.dat")
    dataPath = dataPathHome;
else
    dataPath = dataPathDATA;
end

if kwargs.loadDataFromDATA == 1
    dataPath = P.dataPathDATA + mname + "/" + session_n + "/";
else
    dataPath = P.dataPathHome + mname + "/" + session_n + "/";
end

if kwargs.loadprmFromDATA == 1
    directory = P.pathSpikeSortingDATA;
else
    directory = P.pathSpikeSortingHome;
end

path = fullfile(directory, mname, session_n);

if isfile(path + "Config_h32_oe_bin.prm")
    filename = path + "Config_h32_oe_bin.prm";
elseif isfile(path + "Config_h32_oe_dat.prm")
    filename = path + "Config_h32_oe_dat.prm";
elseif isfile(path + "Config_h32_oe.prm")
    filename = path + "Config_h32_oe.prm";
else
    disp("Params file not found")
    return
end

% Unit construction
directory = P.pathSpikeSortingHome;
config_path = fullfile(directory, mname, P.s(s));
if isfile(fullfile(config_path, "StructEphysData.mat"))
    disp("There is already a unit file");
    continue
end

if ~isfolder(config_path) && s == 2
    disp("Only one session found")
    continue
end

if ~isempty(dir(fullfile(config_path, "*.csv")))
    epdir = dir(fullfile(config_path, "*.csv"));
    filename = fullfile(epdir.folder, epdir.name);
else
    disp("csv file not found")
    continue
end


if kwargs.loadDataFromDATA == 1
    dataPath = fullfile(P.pathEphysDataDATA, mname, P.s(s));
else
    dataPath = fullfile(P.pathEphysDataHome, mname, P.s(s));
end
if ~isfile(fullfile(dataPath,'stimtimes.mat'))
    disp("No stimtimes found, this is a Rhymdata mouse.")
    continue
end

CS_only_trial_idcs = 1:11:length(stimtimes);
trials = 1:length(stimtimes);
CS_only_trial_order = [CS_only_trial_idcs  trials(~ismember(trials,CS_only_trial_idcs))];
unit.stimtimes_trial_order_CS_only_first = CS_only_trial_order;


% computeChannelCandidates
fprintf('Contamination percentage of intervals shorter than 2 ms is: %.4f%%\n', contamination.percentage);

% Parameters
totalSpikes = numel(spikeTimes); % IK change. 18/12/23. made a new way for determining contamination, based on methods by Banga et al., 2023
recordingDuration = max(spikeTimes);
aveFiringRate = totalSpikes/recordingDuration; % Average spiking rate (events per unit time)

maxRefractoryPeriod = 0.002; % ms % Maximum refractory period (time interval)

% Calculate total number of spikes allowed to exist considering the refractory period (for example, within a certain duration)
thresTotalSpikes = aveFiringRate * maxRefractoryPeriod;

% Calculate the probability
prob_threshold = 0.1; % 10% threshold
threshold_spikes = round(prob_threshold * thresTotalSpikes);

% Initialize probability sum
prob = 0;

% Compute cumulative probability using Poisson distribution
for k = 0:threshold_spikes
    %prob = prob + exp(-aveFiringRate * maxRefractoryPeriod) * (aveFiringRate * maxRefractoryPeriod)^k / factorial(k);
    prob = prob + exp(-aveFiringRate * maxRefractoryPeriod) * (aveFiringRate * maxRefractoryPeriod)^k / factorial(k);
end

% Display the probability
fprintf('Probability of observing fewer than %.1f%% of spikes violating the maximum refractory period: %.4f%%\n', prob_threshold * 100, prob*100);

contamination.percentage = 100-prob*100;
%
Parameters

minimum_refractory_period = 0.001; % Minimum refractory period (time interval)

% Calculate adjusted firing rate to respect the minimum refractory period
lambda_adjusted = aveFiringRate * (1 - minimum_refractory_period);

% Compute the probability of observing no spikes within the minimum refractory period
prob_no_spikes = exp(-lambda_adjusted * minimum_refractory_period);

% Calculate the probability of observing fewer than 10% of spikes violating the minimum refractory period
prob_threshold = 0.1; % 10% threshold
threshold_spikes = round(prob_threshold * (1 / lambda_adjusted));

% Initialize cumulative probability
cumulative_prob = prob_no_spikes;

% Compute cumulative probability using Poisson distribution
for k = 1:threshold_spikes
    cumulative_prob = cumulative_prob + exp(-lambda_adjusted) * (lambda_adjusted^k) / factorial(k);
end

% Display the probability
fprintf('Probability of observing fewer than %.1f%% of spikes violating the minimum refractory period: %.4f%%\n', prob_threshold * 100, cumulative_prob*100);

%%
% Calculate total spikes and recording duration
totalSpikes = numel(spikeTimes);
recordingDuration = max(spikeTimes);

% Calculate average firing rate
aveFiringRate = totalSpikes / recordingDuration; % Average spiking rate (events per unit time)

% Maximum refractory period (time interval)
maxRefractoryPeriod = 0.001; % ms

% Calculate the adjusted firing rate considering the maximum refractory period
lambda_adjusted = aveFiringRate * (1 - exp(-aveFiringRate * maxRefractoryPeriod));

% Calculate the probability of observing spikes violating the maximum refractory period
prob_threshold = 0.1; % 10% threshold
threshold_spikes = round(prob_threshold * (1 / lambda_adjusted));

% Initialize cumulative probability
cumulative_prob = 0;

% Compute cumulative probability using Poisson distribution
for k = 0:threshold_spikes
    cumulative_prob = cumulative_prob + exp(-lambda_adjusted) * (lambda_adjusted^k) / factorial(k);
end

% Probability of spikes violating the maximum refractory period
prob_violating_max_refractory = 1 - cumulative_prob;

% Display the probability
fprintf('Probability of observing more than 90%% of spikes violating the maximum refractory period: %.4f%%\n', 100 * prob_violating_max_refractory);

%%
% Calculate total spikes and recording duration
totalSpikes = numel(spikeTimes);
recordingDuration = max(spikeTimes);

% Calculate average firing rate
aveFiringRate = totalSpikes / recordingDuration; % Average spiking rate (events per unit time)

% Minimum refractory period (time interval)
maxRefractoryPeriod = 0.001; % 1 ms

% Calculate adjusted firing rate considering the minimum refractory period
lambda_adjusted = aveFiringRate * (1 - exp(-aveFiringRate * maxRefractoryPeriod));

% Calculate the probability of observing no spikes within the minimum refractory period
prob_no_spikes = exp(-lambda_adjusted * maxRefractoryPeriod);

% Calculate the probability of observing fewer than 10% of spikes violating the minimum refractory period
prob_threshold = 0.1; % 10% threshold
threshold_spikes = round(prob_threshold * (1 / lambda_adjusted));

% Initialize cumulative probability
cumulative_prob = prob_no_spikes;

% Compute cumulative probability using Poisson distribution
for k = 1:threshold_spikes
    cumulative_prob = cumulative_prob + exp(-lambda_adjusted) * (lambda_adjusted^k) / factorial(k);
end

% Probability of observing fewer than 10% of spikes within the minimum refractory period
prob_less_than_10 = cumulative_prob;

% Display the probability
fprintf('Probability of observing fewer than %.1f%% of spikes within the 1 ms refractory period: %.4f%%\n', prob_threshold * 100, 100*prob_less_than_10);



smoothing the lines in the figure
samplingRateIncrease = 100;
newXSamplePoints = linspace(1, size(psthCountsPerCondition, 2), size(psthCountsPerCondition, 2)*samplingRateIncrease);
smoothedData = spline(1:size(psthCountsPerCondition, 2), psthCountsPerCondition, newXSamplePoints);
function smoothed_data = moving_average(data, window_size)
smoothed_data = filter(ones(1, window_size) / window_size, 1, data);
end

% Usage:
% Assuming 'y' contains your data and 'window_size' is the size of the moving average window
smoothedData(1,:) = psthCountsPerCondition(1,:);
smoothedData(2,:) = moving_average(psthCountsPerCondition(2,:), 100);

nConditions = 2; % IK change. TODO: Change later to include US-only trials

P = IkUtils.getParams(); % IK
bins_cs = -(P.t_pre_trial):0.0005:(P.t_post_trial); %IK
rangeMask_new = psthBins >= -0.5 & psthBins <= 2; % IK
iets = psthCountsPerCondition(c,rangeMask(1:end-1));

set(0, 'DefaultFigureRenderer', 'painters');
set(gcf, 'Renderer', 'painters');

% plotRasterWithHistogram

spikeTimes = rasterXY(1,1:3:end) * 1000; % IK change 2 to 3. 8-1-24

binWidth = IkUtils.getParams().BinW;
edges = (range.min:binWidth:range.max) * 1000;

rangeMask = spikeTimes >= range.min * 1000 & spikeTimes <= range.max * 1000;

counts = histcounts(spikeTimes(rangeMask), edges);

maxCount = max(counts);
if any(maxCount == [0, nan, inf])
    normalization = 1;
else
    normalization = maxCount;
end

Highlight regions conditions

lastCSonlyTrial = floor(max(rasterXY_adapted(2, find(rasterXY_adapted(3,:) == 0))));
firstUSonlyTrial = min(rasterXY_adapted(2, find(rasterXY_adapted(3,:) == 2)));
lastUSonlyTrial = floor(max(rasterXY_adapted(2, find(rasterXY_adapted(3,:) == 2))));

% plot summary figures

figure()
boxplot_data = {[allStats.us(non_mod_mask).maxAmpTime], [allStats.cs(cs_facilitation_mask).maxAmpTime], [allStats.cs(cs_suppression_mask).maxAmpTime],[allStats.us(us_facilitation_mask).maxAmpTime], [allStats.us(us_suppression_mask).maxAmpTime]};
boxplot_labels = {'Non modulating','CS facilitation', 'CS suppression', 'US facilitation', 'US suppression'};
hold on
for i = 1:numel(boxplot_data)
    boxplot(boxplot_data{i}, 'positions', i, 'orientation','horizontal')
end
ylabel('Modulation type')
xlabel('Peak Time')
xlim([0 0.35]);
yticks(1:5);
yticklabels(boxplot_labels)
set(gca,'tickdir','out','box','off');
title('Peak time neuron frequency per modulation type')
hold off

% tracesOverview

subplot(2,2,1)
hold on
plot(ts(1:41), squeeze(WTtraces(1,WTtype(1:end)==1,1:41)),'color', [0.5 0.5 0.5])
plot(ts(41:91), squeeze(WTtraces(lastDay,WTtype(1:end)==1,41:91)),'color', [0 0 1])
plot(ts(91:95), squeeze(WTtraces(lastDay,WTtype(1:end)==1,91:95)),'color', [0 1 0])
plot(ts(95:nTimeSteps), squeeze(WTtraces(lastDay,WTtype(1:end)==1,95:nTimeSteps)),'color', [0.5 0.5 0.5])
plot(ts, nanmean(squeeze(WTtraces(lastDay,WTtype(1:end)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: WT', lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

%%
% Source and destination directory paths
% mcode = "MI23.05990.05";
mouseListCopying = ["MI23.03123.02", "MI23.03123.03", "MI23.03123.04", "MI23.03374.01", "MI23.03374.03", "MI23.05990.02", "MI23.05990.03"];
for m = 1 : length(mouseListCopying)
    mcode = mouseListCopying(m);
    mcodeParts = split(mcode, ".");
    if length(mcodeParts) > 1
        part1_ = char(mcodeParts(1));
        part1 = string(part1_(1)) + string(part1_(2));
        NASmcode = part1 + mcodeParts(2) + mcodeParts(3);
    else
        mcodeParts = split(mcode, "-");
        NASmcode = mcodeParts(2) + mcodeParts(3);
    end


    sourceDir = sprintf("/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
    destDir = sprintf('/home/i.klinkhamer/Documents/Data/behaviorData/%s/', mcode);

    secondEphysDay2Done = true;
    secondEphys = 0;

    % List all files in the source directory
    folderList = dir(fullfile(sourceDir)); % You can specify the file extension you want to copy
    folders = "";
    idx = 1;
    for i = 1 : length(folderList)
        if folderList(i).isdir && ~strcmp(folderList(i).name, '.') && ~strcmp(folderList(i).name, '..')
            folders(idx) = string(fullfile(folderList(i).folder, folderList(i).name));
            idx = idx + 1;
        end
    end

    for j = 1:length(folders)
        fileList = dir(folders(j));
        folderExtensions = "";
        for jj = 1 : length(fileList)
            [~, ~, folderExt] = fileparts(fileList(jj).name);
            folderExtensions(jj) = string(folderExt);
        end
        fileList = fileList(~[fileList.isdir] & folderExtensions ~= ".mp4");

        % path_list = IkUtils.getPathFromDir(dir_content);
        % Loop through each file and copy it to the destination directory
        lastSession = 1;
        for ii = 1:numel(fileList)
            parts = split(fileList(ii).name, '_');
            pattern = 's\d+\d+';
            mask = regexp(parts, pattern);
            matching_idx = find(cellfun(@(x) ~isempty(x), mask));
            parts2 = split(parts(matching_idx), "");
            lastSession = max(lastSession, str2double(string(parts2{3}) + string(parts2{4})));
        end

        if j == 12 && secondEphysDay2Done
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession - 1);
                pattern2New = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession - 1);
                pattern2New = sprintf("s%d", lastSession);
            end
            secondEphys = true;
        else
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession);
                pattern2New = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession);
                pattern2New = sprintf("s%d", lastSession);
            end
        end
        partsFile = split([fileList.name], '_');
        pattern = 's\d+\d+';
        maskPattern = regexp(partsFile, pattern);
        matching_indices = find(cellfun(@(x) ~isempty(x), maskPattern));
        newMaskPatternFiles = partsFile(matching_indices);
        newMaskFiles = strcmp(newMaskPatternFiles, patternNew);
        newMaskFiles2 = strcmp(newMaskPatternFiles, pattern2New);
        if sum(newMaskFiles) < 482
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession);
            end
            newMaskFiles = ~strcmp(newMaskPatternFiles, patternNew) & ~newMaskFiles;
        end

        %         maskFiles = strcmp(partsFile, pattern);
        %         maskFiles2 = strcmp(partsFile, pattern2);


        for i = 1:numel(fileList)
            %         if secondEphysDone == true && j == 12
            %             partsFile = split(fileList(i).name, '_');
            %             if lastSession < 10
            %                 pattern = sprintf("s0%d", lastSession - 1);
            %                 pattern2 = sprintf("s0%d", lastSession);
            %             else
            %                 pattern = sprintf("s%d", lastSession - 1);
            %                 pattern2 = sprintf("s%d", lastSession);
            %             end
            %             maskFiles = strcmp(partsFile, pattern);
            %             maskFiles2 = strcmp(partsFile, pattern2);
            %             if sum(maskFiles) < 482
            %                 if lastSession < 10
            %                     pattern = sprintf("s0%d", lastSession);
            %                 else
            %                     pattern = sprintf("s%d", lastSession);
            %                 end
            %                 maskFiles = ~strcmp(partsFile, pattern) & ~maskFiles;
            %             end
            %             secondEphys = 1;
            %         else
            %             partsFile = split(fileList(i).name, '_');
            %             if lastSession < 10
            %                 pattern = sprintf("s0%d", lastSession);
            %             else
            %                 pattern = sprintf("s%d", lastSession);
            %             end
            %             maskFiles = strcmp(partsFile, pattern);
            %         end
            if newMaskFiles(i)%any(maskFiles)
                [~, folderName, ~] = fileparts(folders(j));
                destFolder = fullfile(destDir, folderName);
                if ~exist(destFolder, 'dir')
                    status = mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResults/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    status = mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
            if secondEphys && newMaskFiles2(i)
                [~, folderName, ~] = fileparts(folders(j));
                folderNameNew = sprintf("%s2", folderName);
                destFolder = fullfile(destDir, folderNameNew);
                if ~exist(destFolder, 'dir')
                    status = mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResults/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderNameNew);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    status = mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
        end
        fprintf("folder %d done \n", j)
    end

end
%%
% Specify the directory containing the files
directory = '/home/i.klinkhamer/Documents/Data/spikeSortingUnits/';

% List all files in the directory
% Change '*.txt' to match the file extension of your files

% Specify the new part to insert
newPart = '_20-MI19154-08_2021-05-21'; % Modify this according to your requirements


% Loop through each file and rename it
p = IkUtils.getParams();
mouseList = P.mouseList;

for mname = mouseList
    for s = 1:3
        files = dir(fullfile(directory, mname, p.s(s)));
        if isempty(files)
            continue
        end
        stamp = nameDateStampFiles(mcode = mname, s = s);
        for i = 1:numel(files)
            if ~files(i).isdir
                % Split filename and extension
                [filePath, fileName, fileExt] = fileparts(files(i).name);

                parts = split(fileName, '_');
                if any(strcmp(parts, mname))
                    continue
                end

                % Construct the new filename with the new part inserted between original filename and extension
                newFileName = fullfile(directory, mname, p.s(s), fileName + "_" + stamp + fileExt);

                % Move the file to the new filename
                movefile(fullfile(directory, mname, p.s(s), files(i).name), newFileName);
            end
        end
    end
end

%%
load("/home/i.klinkhamer/PycharmProjects/pythonProject/StructEphysData.mat")
spiketimes = unit.neuron(1).spiketimes;
trials = unit.neuron(1).trial_numbers;


h5create('/home/i.klinkhamer/PycharmProjects/pythonProject/spiketimes.h5', '/dataset1', size(spiketimes));
h5write('/home/i.klinkhamer/PycharmProjects/pythonProject/spiketimes.h5', '/dataset1', spiketimes)
save("/home/i.klinkhamer/PycharmProjects/pythonProject/spiketimes.mat", "spiketimes")
save("/home/i.klinkhamer/PycharmProjects/pythonProject/spiketimes.mat", "spiketimes")
save("/home/i.klinkhamer/PycharmProjects/pythonProject/trials.mat", "trials")
bins = -3:0.0005:3;
trial_index = zeros(length(bins), 220);
trial_index(find(bins >= 0 & bins <= 0.27), :) = 1;

save("/home/i.klinkhamer/PycharmProjects/pythonProject/trials.mat", "trial_index");
%%
trialdata2 = behavior_trial_data;
for i = fieldnames(behavior_trial_data)'
    var = trialdata2.(i{1});
    try
        var(:,:) = nan;
    catch
        [rows, cols] = size(var);
        var = num2cell(NaN(rows, cols));
    end
    trialdata2.(i{1}) = var;
    newfolder = fullfile(fileparts(folder),"230525");
    behavior_trial_data_save = behavior_trial_data;
    behavior_trial_data = trialdata2;
    save(fullfile(newfolder, 'trialdata.mat'), 'behavior_trial_data');
    behavior_trial_data = behavior_trial_data_save;
end

%%
fprintf("\n-------------------------------------")
fprintf("\nWhich mouse would you like to curate?")
fprintf("\n-------------------------------------\n")

mouseCodes = arrayfun ...          % IK change
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
mouseNames = arrayfun ...
    ( @(mouse) string(mouse.name) ...
    , defaultMice() ...
    );
%mouseNames = [mouseNames, "All of the above"];


[mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
mcode = mouseCodes(mname_idx);
%mname = mouseNames(mname_idx);

p = IkUtils.getParams();

%path_trial_data = which(fullfile(mcode, 'trialdata.mat'));
dir_content = dir(fullfile(p.pathBehaviorDataHome, mcode));
path_list = IkUtils.getPathFromDir(dir_content);
% dates = dir_content(3:end).name;

for f = 1 : length(path_list)

    folder = path_list(f);
    try
        load(fullfile(folder, 'trialdata.mat'));
    catch
        continue
    end
    continuefromnow = false;
    lastTrial = length(behavior_trial_data.eyelidpos);

    if lastTrial < 240
        trialdata2 = behavior_trial_data;
        for i = fieldnames(behavior_trial_data)'
            var = trialdata2.(i{1});
            new_var = var;
            [rows, cols] = size(var);
            if i == "meanAll"
                continuefromnow = true;
            end
            lastTrial = length(behavior_trial_data.eyelidpos);
            if i == "CRamp"
                continuefromnow = false;
                lastTrial = length(trialdata2.CRamp);
            end
            if continuefromnow
                continue
            end
            if isnumeric(var) && isequal(class(var), 'double') && i ~= "CRamp5"
                if rows == 1 || rows == 200
                    new_var = nan(rows, 240);
                    new_var(:,1:lastTrial) = var;
                else
                    new_var = nan(240, cols);
                    if i == "CRamp"
                        new_var = nan(220,cols);
                    end
                    new_var(1:lastTrial,:) = var;
                end
                %                 try
                %
                %                     new_var(lastTrial,:) = nan;
                %                 catch
                %                     new_var(:,lastTrial) = nan;
                %                 end
            elseif isequal(class(var), 'string')
                new_var = var;
                new_var(lastTrial+1:240) = '';
            elseif isequal(class(var), 'cell')
                [rows, cols] = size(var);
                %                 new_var = cell(rows, cols);
                %new_var{:,:} = nan;
                if rows == 1 || rows == 200
                    new_var = cell(rows, 240);
                    for ii = 1:rows
                        for j = 1:lastTrial
                            new_var{ii, j} = var{ii,j};
                        end
                    end
                elseif cols == 1 || cols == 200
                    new_var = cell(240, cols);
                    for ii = 1:lastTrial
                        for j = 1:cols
                            new_var{ii, j} = var{ii, j};
                        end
                    end
                end

            end
            trialdata2.(i{1}) = new_var;
            %newfolder = fullfile(fileparts(folder),"230525");
            if i == "CRamp"
                continuefromnow = true;
            end
        end
        save(fullfile(folder, 'trialdataor.mat'), 'behavior_trial_data');
        behavior_trial_data_save = behavior_trial_data;
        behavior_trial_data = trialdata2;
        save(fullfile(folder, 'trialdata.mat'), 'behavior_trial_data');
        behavior_trial_data = behavior_trial_data_save;
        continuefromnow = false;
    end
end
%%
if isfile(fullfile(path, '100_CH1.continuous'))
    formatSpec = '%s%d_CH%d.continuous'; %filename format for .continuous data
    format = '100_CH';
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
elseif isfile(fullfile(path, '100_1.continuous'))
    formatSpec = '%s%d_%d.continuous'; %filename format for .continuous data
    format = '100_';
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
elseif isfile(fullfile(path, '100_CH1_2.continuous'))
    formatSpec = '%s%d_CH%d_2.continuous'; %filename format for .continuous data
    disp("Format contains 2 and I haven't fixed a way to deal with that in the code yet, so this is going to create an error.")
    format = '';
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
else
    disp("No channel.continuous file found")
    filepath = [];
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
    channel = chSel(ch);
    filepath = fullfile(path, sprintf('100_CH%s.continuous', string(channel)));
    %filepath = fullfile(path, '100', (channel), '.continuous');
    return
end
%%
function reSaveWithStamp(filecode, filetype, mname, s)
arguments
    filecode = "StructEphysData"
    filetype = ".mat"
    mname = "20-MI19442-05"
    s = 1
end

P = IkUtils.getParams();
mouseList = P.mouseList;

for mname = mouseList
    for s = 1:3

        fprintf("Mouse: %s, %s  Session: %d \n", mname, P.mouseNames(find(P.mouseList == mname)), s);

        [filename, folder] = recursiveFileSearch(filecode + "*" + filetype, mname, s);
        if length(filename) ~= 1
            continue
        end

        structFilename = filecode + filetype;

        fullfilename = fullfile(folder, filename);
        path = fileparts(fullfilename);

        stamp = nameDateStampFiles(mcode = mname, s = s, file_pattern = filecode, file_extension = filetype);
        loaded_data = load(fullfilename);
        var_name = fieldnames(loaded_data);
        var_name = var_name{1};
        new_variable = loaded_data.(var_name);
        unit = loaded_data.(var_name);
        %         load(fullfilename)

        save(fullfile(path, stamp),'unit')
    end
end




% filename_pattern = sprintf('StructEphysData_*\d{8}.mat'); % filename_pattern = sprintf('StructEphysData_*.mat');
% matching_files = dir(fullfile(path, filename_pattern));
% % Display the list of matching files
% disp({matching_files.name});


end

%%
filepath = fullfile(path, sprintf('100_%s%d%s%s.continuous', formatPart, ch, formatExtra));
if parts{2} == "CH1"
    if parts{3} == "2"
        format = '100_CH%d_2';
    else
        format = '100_CH%d';
    end
    if ch > 16
        ch = ch+32;
    end
elseif parts{2} == "1"
    format = '100_';
else
    disp("No channel.continuous file found")
    return
end
filepath = fullfile(path, sprintf('%s%d.continuous', format, ch));

%%
function reSaveWithStamp(kwargs)
arguments
    kwargs.directory = '/home/i.klinkhamer/Documents/Data/spikeSortingUnits/';
    kwargs.AllSessions = true;
    kwargs.AllMice = true;
    kwargs.startMouse = 1;
end
directory = kwargs.directory;
P = IkUtils.getParams();
if kwargs.AllMice == 0

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse's data would you like resave?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );

    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);
    mousenames = mcode;
else
    mousenames = P.mouseList(kwargs.startMouse:end);
end

if kwargs.AllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    sessions = 1:P.n_sessions;
end

% Loop through each file and rename it
for mname = mousenames
    %     prevDate = "";
    secondEphysSessionOfDay2 = 0;
    for s = sessions
        files = dir(fullfile(directory, mname, P.s(s)));
        if isempty(files)
            continue
        end
        stamp = nameDateStampFiles(mcode = mname, s = s);
        %         currDate = split(stamp, '_'); currDate = currDate{2};
        %         try
        %             if currDate == prevDate
        %                 stamp = stamp + "_2";
        %                 secondEphysSessionOfDay2 = 1;
        %             end
        %         catch
        %         end
        for i = 1:numel(files)
            if ~files(i).isdir
                % Split filename and extension
                [filePath, fileName, fileExt] = fileparts(files(i).name);

                parts = split(fileName, '_');
                if any(strcmp(parts, mname)) && ~secondEphysSessionOfDay2
                    continue
                    %                 elseif any(strcmp(parts, mname)) && any(strcmp(parts, "2")) && secondEphysSessionOfDay2
                    %                     continue
                    %                 elseif any(strcmp(parts, mname)) && ~any(strcmp(parts, "2")) && secondEphysSessionOfDay2
                    %                     stamp = "2";
                end


                % Construct the new filename with the new part inserted between original filename and extension
                newFileName = fullfile(directory, mname, P.s(s), fileName + "_" + stamp + fileExt);

                % Move the file to the new filename
                movefile(fullfile(directory, mname, P.s(s), files(i).name), newFileName);
            end
        end
        %         prevDate = currDate;
        %         secondEphysSessionOfDay2 = 0;
    end
end

end

%%
%% IK 2-11-23
function batchProcessTrialsExampleIK(mInput)
close all
debug = 0;
seeFigures = 1;
plotfigures = 0;
% trainingSessions = true;
p = IkUtils.getParams();
%addPaths()
% addpath("/home/mick/Desktop/Ilse/ephys_code/");
% addpath("/home/mick/Desktop/Ilse/ephys_code/Complex spike suite revised/");
% addpath('/media/mick/DATA/Ilse/Data/convertedData/')
% addpath("/home/mick/Desktop/Ilse/ephys_code/Code Nynke/analysis code/")
% addpath("/home/mick/Desktop/Ilse/ephys_code/Code Nynke/neuroblinks-master/utilities/")
%base_dir = 'C:\Users\nmhet\Documents\studie\MEP\data';
mouseList = p.mouseList;



if nargin < 1
    %     mcode = "20-MI19442-03";

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to curate?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...          % IK change
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );
    %mouseNames = [mouseNames, "All of the above"];


    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);
    %mname = mouseNames(mname_idx);
    % mcode = "MI0599008";
else
    if any(contains(mInput, p.mouseNames))
        mcode = p.mouseList(find(p.mouseNames == mInput));
    else
        mcode = mInput;
    end
end



% pathNAS = "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/";
% pathNAS = "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/";
% path_trial_data = which(fullfile(mcode, 'trialdata.mat'));

% dir_content = dir(fullfile(pathNAS, mcode));
try
    dir_content = dir(fullfile(p.pathBehaviorDataHome, mcode));
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    trainingPath = "/home/i.klinkhamer/Documents/Data/behaviorDataTrainingOnlyMice/";
    dir_content = dir(fullfile(trainingPath, mcode));
    path_list = IkUtils.getPathFromDir(dir_content);
    trainingOnly = true;
end
% dates = dir_content(3:end).name;

for f = 1 : length(path_list)

    folder = path_list(f);
    [folderPath, folderName, ~] = fileparts(folder);
    %             folder = fullfile(p.pathBehaviorDataHome, mcode, folderName);
    if ~exist(folder, 'dir')
        status = mkdir(folder);
    end
    if ~trainingOnly
        behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResults/";
    else
        behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
    end
    behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
    if ~exist(behaviorDataResultsFolder, 'dir')
        status = mkdir(behaviorDataResultsFolder);
    end
    s = find(p.s == erase(erase(folder, fileparts(folder)),"/"));
    if isempty(s)
        s = erase(erase(folder, fileparts(folder)),"/");
    end
    stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyelidpos', file_extension = '.mat');
    stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');

    if (~exist(fullfile(folder, 'trialdata.mat'), 'file') && ~exist(fullfile(folder, stamptrial), 'file')) || debug == 1
        behavior_trial_data = processTrials(folder, 'recalibrate'); %,'recalibrate');  % Recalibrate eyelid
    elseif seeFigures == 1
        try
            load(fullfile(folder, stamptrial));
        catch
            load(fullfile(folder, 'trialdata.mat'));
        end
        %behavior_trial_data = behavior_trial_data.behavior_trial_data;
    else
        continue
    end
    %             try
    %                 behavior_trial_data = [behavior_trial_data.behavior_trial_data];
    %                 try
    %                     behavior_trial_data = [behavior_trial_data.behavior_trial_data];
    %                     save(fullfile(folder, 'trialdata.mat'), 'behavior_trial_data');
    %                 catch
    %                     save(fullfile(folder, 'trialdata.mat'), 'behavior_trial_data');
    %                 end
    %             catch
    %             end

    if ~isempty(behavior_trial_data)

        eyelidpos = behavior_trial_data.eyelidpos;
        min_val_eye = min(eyelidpos(:));
        max_val_eye = max(eyelidpos(:));
        eyelidpos = (eyelidpos - min_val_eye) / (max_val_eye - min_val_eye);
        %         eyelidpos_upsampled = repelem(eyelidpos, 1,10);
        % Define the factor by which you want to upsample
        upsample_factor = 10;

        % Create the original time points
        original_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2));

        % Create the new time points after upsampling
        new_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2)*upsample_factor);

        % Interpolate the data
        eyelidpos_upsampled = interp1(original_time, eyelidpos.', new_time, 'spline');
        % Gaussian smoothing
        %         sigma = 5; % Standard deviation for Gaussian smoothing
        %         eyelidpos_smoothed = imgaussfilt(eyelidpos_upsampled, sigma);

        %         eyelidpos = eyelidpos.';
        %         eyelidpos_smoothed = eyelidpos_smoothed.';

        min_val_eye = min(eyelidpos_upsampled(:));
        max_val_eye = max(eyelidpos_upsampled(:));
        eyelidpos = (eyelidpos_upsampled - min_val_eye) / (max_val_eye - min_val_eye);
        %         figure;
        %         times = -0.2:0.0005:0.7995;
        %         plot(times, eyelidpos)
        % figure;
        % times = -0.2:0.0005:0.7995;
        % plot(times, eyelidpos)


        save(fullfile(folder, stampeye), 'eyelidpos');
        if (~exist(fullfile(folder, 'trialdata.mat'), 'file') || ~exist(fullfile(folder, stamptrial), 'file')) || debug == 1
            save(fullfile(folder, stamptrial), 'behavior_trial_data');
        end

        try
            first_part_folder = fileparts(fileparts(fileparts(folder)));
            middle_folder = 'behaviorDataResults';
            end_folder = erase(folder, fileparts(fileparts(folder)));
            new_folder = fullfile(first_part_folder,middle_folder,end_folder);
            save(fullfile(new_folder, stampeye), 'eyelidpos');

            if ~exist(fullfile(new_folder, stamptrial), 'file') || debug == 1
                save(fullfile(new_folder, stamptrial), 'behavior_trial_data');
            end
        catch
            disp("not saved in behavior results folder")
        end

        baseline = mean(mean(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,1:40)));

        fullclosure = mean(max(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,40:95),[],2));
        traces = (behavior_trial_data.eyelidpos - baseline) ./ (fullclosure - baseline);
        %         min_val = min(traces(:));
        %         max_val = max(traces(:));
        %         traces = (traces - min_val) / (max_val - min_val);

        %         Assuming behavior_trial_data is a table or array with eyelidpos column
        %         figure
        %         plot(behavior_trial_data.tm(1,:), traces)
        %         title("original")
        %
        %         figure
        %         plot(behavior_trial_data.tm(1,:), behavior_trial_data.tracesnorm)
        %         title("original norm")

        traces_new = behavior_trial_data.eyelidpos - (mean(behavior_trial_data.eyelidpos(:,1:40), 2) - baseline);
        baseline_new = mean(mean(traces_new(behavior_trial_data.c_csdur==0,1:40)));
        stdev = std(behavior_trial_data.eyelidpos(:,1:40), 0, 2);
        mask_stdev = stdev < 3*mean(stdev);
        mask_lowest = min(traces_new(:,41:end), [], 2) > (baseline_new - 0.3);
        mask_highest = max(traces_new(:,141:end), [],2) < max(max(traces_new(:,90:140)));
        traces_new = traces_new(mask_stdev & mask_lowest & mask_highest, :);
        min_val = min(traces_new(:));
        max_val = max(traces_new(:));
        traces_new = (traces_new - min_val) / (max_val - min_val);

        if plotfigures
            figure; hold on
            x1 = 0;
            x2 = 250;
            y1 = 0;
            y2 = 1;
            patch ...
                ( ...
                [x1 x1 x2 x2] ...
                , [y1 y2 y2 y1] ...
                , 'blue' ...
                , EdgeColor = 'none' ...
                , FaceAlpha = 0.2 ...
                )
            x1 = 250;
            x2 = 270;
            patch ...
                ( ...
                [x1 x1 x2 x2] ...
                , [y1 y2 y2 y1] ...
                , 'green' ...
                , EdgeColor = 'none' ...
                , FaceAlpha = 0.2 ...
                )
            plot(behavior_trial_data.tm(1,:), traces_new)

            title("Eyelid closure")
        end
    end
end

end
%%
for i=1:size(WTmice,2)
    mouse = WTmice(i);
    dir_content = dir(fullfile(base_dir, mouse));
    path_list = IkUtils.getPathFromDir(dir_content);
    for f = 1 : length(path_list)
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mouse, s, folder);
        load(fullfile(folder, filename(end)));
        trial_data(f) = behavior_trial_data;
        traces_new(f).traces = normalizeEyelidTraces(behavior_trial_data);
        trial_data(f).tracesnorm = traces_new(f).traces;
        %trial_data(f).traces_new = traces_new;
    end
    trial_data_old = trial_data;

    trial_data = trial_data_recalc(trial_data);

    if length(trial_data(f).CRamp)<220
        trial_data(f).CRamp(length(trial_data(f).CRamp)+1:220) = nan;
    end



    for j=1:p.lastDay %length(trial_data)
        WTtraces(j,(i*p.nTrials-(p.nTrials-1)):i*p.nTrials,:) = trial_data(j).tracesnorm(1:p.nTrials,:);
        WTcramp(j,i,1:min(nCSplusRegTrials, length(trial_data(j).CRamp))) = trial_data(j).CRamp(1:min(nCSplusRegTrials, length(trial_data(j).CRamp)))';
        WTcramp5(j,i,1:length(trial_data(j).CRamp5)) = trial_data(j).CRamp5';
        WTCRperc(j,i) = sum(~isnan(WTcramp5(j,i,:)))/nCSplusRegTrials*100;

        WT(i).s(j).traces = trial_data(j).tracesnorm(1:p.nTrials,:);
        WT(i).s(j).cramp = trial_data(j).CRamp(1:min(nCSplusRegTrials, length(trial_data(j).CRamp)))';
        WT(i).s(j).cramp5 = trial_data(j).CRamp5';
        WT(i).s(j).CRperc = sum(~isnan(WTcramp5(j,i,:)))/nCSplusRegTrials*100;
    end
end
WTCRperc(WTCRperc == 0) = nan;


for i=1:size(MUTmice,2)
    mouse = MUTmice(i);
    dir_content = dir(fullfile(base_dir, mouse));
    path_list = IkUtils.getPathFromDir(dir_content);
    for f = 1 : length(path_list)
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mouse, s, folder);
        load(fullfile(folder, filename(end)));
        %         load(fullfile(folder, 'trialdata.mat'));
        trial_data(f) = behavior_trial_data;
        traces_new(f).traces = normalizeEyelidTraces(behavior_trial_data);
        %trial_data(f).traces_new = traces_new;
        trial_data(f).tracesnorm = traces_new(f).traces;
    end
    trial_data_old = trial_data;

    trial_data = trial_data_recalc(trial_data);

    if length(trial_data(f).CRamp)<220
        trial_data(f).CRamp(length(trial_data(f).CRamp)+1:220) = nan;
    end
    for j=1:p.lastDay%length(trial_data)
        MUTtraces(j,(i*p.nTrials-(p.nTrials-1)):i*p.nTrials,:) = trial_data(j).tracesnorm(1:p.nTrials,:);
        MUTcramp(j,i,1:min(nCSplusRegTrials, length(trial_data(j).CRamp))) = trial_data(j).CRamp(1:min(nCSplusRegTrials, length(trial_data(j).CRamp)))';
        MUTcramp5(j,i,1:length(trial_data(j).CRamp5)) = trial_data(j).CRamp5';
        MUTCRperc(j,i) = sum(~isnan(MUTcramp5(j,i,:)))/nCSplusRegTrials*100;

        MUT(i).s(j).traces = trial_data(j).tracesnorm(1:p.nTrials,:);
        MUT(i).s(j).cramp = trial_data(j).CRamp(1:min(nCSplusRegTrials, length(trial_data(j).CRamp)))';
        MUT(i).s(j).cramp5 = trial_data(j).CRamp5';
        MUT(i).s(j).CRperc = sum(~isnan(MUTcramp5(j,i,:)))/nCSplusRegTrials*100;
    end
end
MUTCRperc(MUTCRperc == 0) = nan;

type(trial_data(1).c_csdur == 0) = 0;
type(trial_data(1).c_usdur == 0) = 2;
WTtype = zeros(length(WTmice), p.nTrials);
MUTtype = zeros(length(MUTmice), p.nTrials);
for i = 1 : length(WTmice)
    WTtype(i, :) = type;
    WT(i).type = type;
end
for i = 1 : length(MUTmice)
    MUTtype(i,:) = type;
    MUT(i).type = type;
end


for j =1:p.lastDay
    WTcrampmean(j)=nanmean(WTcramp(j,:,:), 'all');
    WTcramp5mean(j)=nanmean(WTcramp5(j,:,:), 'all');
    WTcrampstd(j) = nanstd(WTcramp(j,:,:),0, 'all');
    WTcramp5std(j) = nanstd(WTcramp5(j,:,:),0,'all');
    WTmean(j,:) = nanmean(WTtraces(j,WTtype(1:end) == 1,:),2);

    MUTcrampmean(j)=nanmean(MUTcramp(j,:,:), 'all');
    MUTcramp5mean(j)=nanmean(MUTcramp5(j,:,:), 'all');
    MUTcrampstd(j) = nanstd(MUTcramp(j,:,:),0, 'all');
    MUTcramp5std(j) = nanstd(MUTcramp5(j,:,:),0,'all');
    MUTmean(j,:) = nanmean(MUTtraces(j,MUTtype(1:end) == 1,:),2);

    %     WTcrampmean(j)=nanmean(WTcramp(j,:,:), 'all');
    %     WTcramp5mean(j)=nanmean(WTcramp5(j,:,:), 'all');
    %     WTcrampstd(j) = nanstd(WTcramp(j,:,:),0, 'all');
    %     WTcramp5std(j) = nanstd(WTcramp5(j,:,:),0,'all');
    %     WTmean(j,:) = nanmean(WTtraces(j,WTtype(j,:) == 1,:),2);
    %
    %     MUTcrampmean(j)=nanmean(MUTcramp(j,:,:), 'all');
    %     MUTcramp5mean(j)=nanmean(MUTcramp5(j,:,:), 'all');
    %     MUTcrampstd(j) = nanstd(MUTcramp(j,:,:),0, 'all');
    %     MUTcramp5std(j) = nanstd(MUTcramp5(j,:,:),0,'all');
    %     MUTmean(j,:) = nanmean(MUTtraces(j,MUTtype(j,:) == 1,:),2);
end

WTCRperdstd = std(WTCRperc,1,2);
MUTCRperdstd = std(MUTCRperc,1,2);
ts=trial_data(1).tm(1,:);

%CR onset
for i = 1:p.lastDay
    for j=1:(length(WTtraces))
        if ~isnan(WTtraces(i,j,1))
            baseline = nanmean(WTtraces(i,j,1:40));
            stdbase = nanstd(WTtraces(i,j,1:40));

            if find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1)
                WTonsets(i,j) = ts(find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1));
            end
        end
    end
end

WTonsets = WTonsets(:,WTtype(1:end) ==1 |WTtype(1:end) ==2);
%WTonsets(WTcramp <0.05) = NaN;

for i = 1:p.lastDay
    for j=1:(length(MUTtraces))
        %if ~isnan(MUTtraces(i,j,1))
        baseline = nanmean(MUTtraces(i,j,1:40));
        stdbase = nanstd(MUTtraces(i,j,1:40));

        if find(MUTtraces(i,j,:)>(baseline+10*stdbase) , 1)>0
            MUTonsets(i,j) = ts(find(MUTtraces(i,j,:)>(baseline+10*stdbase) , 1));
        end
        %end
    end
end

MUTonsets = MUTonsets(:,MUTtype(1:end) ==1 |MUTtype(1:end) ==2);

MUTonsets(MUTcramp <0.05) = NaN;
%% plotting
meanTraces = figure('Color', 'white', 'Position', [0 0 1500 2000]);
subplot(1,2,1)
plot(ts, WTmean)
legend('day 1', 'day 2','day 3','day 4','day 5','day 6','day 7','day 8','day 9','day 10')
axis([ts(1) ts(end) min(min(WTmean)) max(max(WTmean))])
title('Mean eyelid traces - WT');
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')
subplot(1,2,2)
plot(ts, MUTmean)
legend('day 1', 'day 2','day 3','day 4','day 5','day 6','day 7','day 8','day 9','day 10')
axis([ts(1) ts(end) min(min(MUTmean)) max(max(MUTmean))])
title('Mean eyelid traces - MUT');
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

CRamp = figure('Color', 'white');
hold on
errorbar(WTcramp5mean,WTcramp5std, 'b');
errorbar(MUTcramp5mean,MUTcramp5std, 'r');
xticks(1:10);
ylim([-0.2 1]);
a = get(gca,'XTickLabel');
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR amplitude')
title('CS-US & CS-only trials - CR>5%')
box('off');
legend('WT', 'MUT')

CRamp5 = figure('Color', 'white');
hold on
errorbar(WTcrampmean,WTcrampstd, 'b');
errorbar(MUTcrampmean,MUTcrampstd, 'r');
xticks(1:10);
ylim([-0.2 1]);
a = get(gca,'XTickLabel');
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR amplitude')
title('CS-US & CS-only trials ')
box('off');
legend('WT', 'MUT')

type = ones(1,p.nTrials);
WTtraces=nan(p.lastDay,size(WTmice,2)*p.nTrials,p.nTimeSteps);
% WTtype=ones(p.lastDay,size(WTmice,1)*p.nTrials);
WTcramp=nan(p.lastDay,size(WTmice,2),nCSplusRegTrials);
WTcramp5=nan(p.lastDay,size(WTmice,2),nCSplusRegTrials);
WTmean=nan(p.lastDay,p.nTimeSteps);
WTonsets=nan(p.lastDay,size(WTmice,2)*p.nTrials);
MUTtraces=nan(p.lastDay,size(MUTmice,2)*p.nTrials,p.nTimeSteps);
% MUTtype=ones(p.lastDay,size(MUTmice,1)*p.nTrials);
MUTcramp=nan(p.lastDay,size(MUTmice,2),nCSplusRegTrials);
MUTcramp5=nan(p.lastDay,size(MUTmice,2),nCSplusRegTrials);
MUTmean=nan(p.lastDay,p.nTimeSteps);
MUTonsets=nan(p.lastDay,size(MUTmice,2)*p.nTrials);

%%
%% IK 15-3-2024
clear all
close all
%% preallocation
p = IkUtils.getParams();
loc = "KO";
twoGroups = false;

MUTmice = getCodesMouseGroup("Tsc1KOMUT");
WTL7mice = getCodesMouseGroup("Tsc1L7WT");
MUTL7mice = getCodesMouseGroup("Tsc1L7MUT");
WTmice = getCodesMouseGroup("Tsc1KOWT");


if twoGroups
    if loc == "KO"
        WTmice = WTmice;
        MUTmice = MUTmice;
    elseif loc == "L7"
        WTmice = WTL7mice
        MUTmice = MUTL7mice
    end
else

end

%% calculations
WT = calcMice(WTmice, p);
MUT = calcMice(MUTmice, p);
WTL7 = calcMice(WTL7mice, p);
MUTL7 = calcMice(MUTL7mice, p);


%% Final figures


if twoGroups
    eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);
    reg_trials_mask = repmat(WT.type(1,1:p.nTrials)==1, 1, length(WTmice));
    subplot(2,2,1)
    hold on
    plot(ts(1:41), squeeze(WT.traces(1,WT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(ts(41:91), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(ts(91:95), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(ts(95:p.nTimeSteps), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(ts, nanmean(squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: WT', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')

    subplot(2,2,2)
    hold on
    MUTtracesPlot = MUT.traces(:,:,:);%MUT.traces(:,481:720,:);
    plot(ts(1:41), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(ts(41:91), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(ts(91:95), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(ts(95:p.nTimeSteps), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(ts, nanmean(squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: MUT', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')

    subplot(2,2,3)
    WTmeanWF=WT.mean;
    WTmeanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
    WTmeanWF(:,p.nTimeSteps)=NaN;
    colormap(winter)
    for i=1:p.lastDay
        color(i,1:p.nTimeSteps)=i*(length(WTmeanWF)); %give each line a different color (instead of gradient for heigth profile)
    end
    h1=waterfall(ts,1:p.lastDay,WTmeanWF,color)
    h1.LineWidth=2;
    h1.FaceAlpha=0; %transparent faces
    view([-0.5 -1.5 1]) %determines view angle
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    xlim([-200 800]);
    title(sprintf('Average traces day 1-%d: WT', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    subplot(2,2,4)
    MUTmeanWF=MUT.mean;
    MUTmeanWF(:,1)=NaN;
    MUTmeanWF(:,p.nTimeSteps)=NaN;
    h2=waterfall(ts,1:p.lastDay,MUTmeanWF,color)
    h2.LineWidth=2;
    h2.FaceAlpha=0;
    view([-0.5 -1.5 1])
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    title(sprintf('Average traces day 1-%d: MUT', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

    subplot(2,2,1)
    hold on
    errorbar(WT.crampmean,WT.crampstd, 'b');
    errorbar(MUT.crampmean,MUT.crampstd, 'r');
    xticks(1:p.lastDay);
    ylim([-0.2 1]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR amplitude')
    title('CS-US & CS-only trials')
    box('off');
    legend('WT', 'MUT')

    subplot(2,2,2)
    hold on
    errorbar(nanmean(WT.CRperc,2),nanstd(WT.CRperc,1,2), 'b');
    errorbar(nanmean(MUT.CRperc,2),nanstd(MUT.CRperc,1,2), 'r');
    xticks(1:p.lastDay);
    ylim([0 100]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR percentage')
    title('CS-US & CS-only trials')
    box('off');

    subplot(2,2,3)
    binEdges = 0:25:350;
    %CR onset time
    IkUtils.histogramIK(WT.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(WT.onsets(p.lastDay,:), binEdges) %IkUtils.histogramIK(WT.onsets(p.lastDay-1,:))
    xlim([0 350])
    title('WT onset latency')
    legend('day 1', sprintf('day%d', p.lastDay))

    subplot(2,2,4)
    %heatmap/density plot
    IkUtils.histogramIK(MUT.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(MUT.onsets(p.lastDay,:),binEdges) %IkUtils.histogramIK(MUT.onsets(p.lastDay-1,:)) % right now day 9 because better results...
    xlim([0 350])
    title('MUT onset latency')
else
    eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);
    reg_trials_mask = repmat(WT.type(1,1:p.nTrials)==1, 1, length(WTmice));
    subplot(2,2,1)
    hold on
    plot(WT.ts(1:41), squeeze(WT.traces(1,WT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(WT.ts(41:91), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(WT.ts(91:95), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(WT.ts(95:p.nTimeSteps), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(WT.ts, nanmean(squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: WT', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')

    subplot(2,2,2)
    hold on
    MUTtracesPlot = MUT.traces(:,:,:);%MUTtraces(:,481:720,:);
    plot(WT.ts(1:41), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(WT.ts(41:91), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(WT.ts(91:95), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(WT.ts(95:p.nTimeSteps), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(WT.ts, nanmean(squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: MUT', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')


    subplot(2,2,3)
    WTmeanWF=WT.mean;
    WTmeanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
    WTmeanWF(:,p.nTimeSteps)=NaN;
    colormap(winter)
    for i=1:p.lastDay
        color(i,1:p.nTimeSteps)=i*(length(WTmeanWF)); %give each line a different color (instead of gradient for heigth profile)
    end
    h1=waterfall(WT.ts,1:p.lastDay,WTmeanWF,color);
    h1.LineWidth=2;
    h1.FaceAlpha=0; %transparent faces
    view([-0.5 -1.5 1]) %determines view angle
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    xlim([-200 800]);
    title(sprintf('Average traces day 1-%d: WT', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    subplot(2,2,4)
    MUTmeanWF=MUT.mean;
    MUTmeanWF(:,1)=NaN;
    MUTmeanWF(:,p.nTimeSteps)=NaN;
    h2=waterfall(WT.ts,1:p.lastDay,MUTmeanWF,color);
    h2.LineWidth=2;
    h2.FaceAlpha=0;
    view([-0.5 -1.5 1])
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    title(sprintf('Average traces day 1-%d: MUT', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    fname = 'eyeblinkTracesOverviewTsc1KO.eps';
    saveFigures(eyeblinkTracesOverview, fname);
    fname2 = 'eyeblinkTracesOverviewTsc1KO.png';
    saveFigures(eyeblinkTracesOverview, fname2);



    eyeblinkTracesOverviewL7 = figure('Color', 'white','Position', [10 0 1000 700]);
    subplot(2,2,1)
    hold on
    WTL7tracesPlot = WTL7.traces(:,:,:);%WTL7.traces(:,481:720,:);
    plot(WT.ts(1:41), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(WT.ts(41:91), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(WT.ts(91:95), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(WT.ts(95:p.nTimeSteps), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(WT.ts, nanmean(squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: WTL7', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')

    subplot(2,2,2)
    hold on
    MUTL7tracesPlot = MUTL7.traces(:,:,:);%MUTL7.traces(:,481:720,:);
    plot(WT.ts(1:41), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
    plot(WT.ts(41:91), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
    plot(WT.ts(91:95), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
    plot(WT.ts(95:p.nTimeSteps), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
    plot(WT.ts, nanmean(squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
    title(sprintf('All traces day %d: MUTL7', p.lastDay))
    xlabel('Time from CS (ms)')
    ylabel('fraction eyelid closure')
    subplot(2,2,3)
    WTL7meanWF=WTL7.mean;
    WTL7meanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
    WTL7meanWF(:,p.nTimeSteps)=NaN;
    colormap(winter)
    for i=1:p.lastDay
        color(i,1:p.nTimeSteps)=i*(length(WTL7meanWF)); %give each line a different color (instead of gradient for heigth profile)
    end
    h1=waterfall(WT.ts,1:p.lastDay,WTL7meanWF,color);
    h1.LineWidth=2;
    h1.FaceAlpha=0; %transparent faces
    view([-0.5 -1.5 1]) %determines view angle
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    xlim([-200 800]);
    title(sprintf('Average traces day 1-%d: WTL7', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    subplot(2,2,4)
    MUTL7meanWF=MUTL7.mean;
    MUTL7meanWF(:,1)=NaN;
    MUTL7meanWF(:,p.nTimeSteps)=NaN;
    h2=waterfall(WT.ts,1:p.lastDay,MUTL7meanWF,color);
    h2.LineWidth=2;
    h2.FaceAlpha=0;
    view([-0.5 -1.5 1])
    zlim=([0 1.2]);
    yticks(1:p.lastDay);
    title(sprintf('Average traces day 1-%d: MUTL7', p.lastDay))
    xlabel('Time from CS (ms)')
    zlabel('fraction eyelid closure')
    ylabel('training day')

    fname = 'eyeblinkTracesOverviewTsc1L7.eps';
    saveFigures(eyeblinkTracesOverviewL7, fname);
    fname2 = 'eyeblinkTracesOverviewTsc1L7.png';
    saveFigures(eyeblinkTracesOverviewL7, fname2);

    CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

    subplot(2,2,1)
    hold on
    errorbar(WT.crampmean,WT.crampstd, 'b');
    errorbar(MUT.crampmean,MUT.crampstd, 'r');
    xticks(1:p.lastDay);
    ylim([-0.2 1]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR amplitude')
    title('CS-US & CS-only trials')
    box('off');
    legend('WT KO', 'MUT KO')

    subplot(2,2,2)
    hold on
    errorbar(nanmean(WT.CRperc,2),nanstd(WT.CRperc,1,2), 'b');
    errorbar(nanmean(MUT.CRperc,2),nanstd(MUT.CRperc,1,2), 'r');
    xticks(1:p.lastDay);
    ylim([0 100]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR percentage')
    title('CS-US & CS-only trials')
    box('off');


    subplot(2,2,3)
    hold on
    errorbar(WTL7.crampmean,WTL7.crampstd, 'b');
    errorbar(MUTL7.crampmean,MUTL7.crampstd, 'r');
    xticks(1:p.lastDay);
    ylim([-0.2 1]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR amplitude')
    title('CS-US & CS-only trials')
    box('off');
    legend('WT L7', 'MUT L7')

    subplot(2,2,4)
    hold on
    errorbar(nanmean(WTL7.CRperc,2),nanstd(WTL7.CRperc,1,2), 'b');
    errorbar(nanmean(MUTL7.CRperc,2),nanstd(MUTL7.CRperc,1,2), 'r');
    xticks(1:p.lastDay);
    ylim([0 100]);
    a = get(gca,'XTickLabel');
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR percentage')
    title('CS-US & CS-only trials')
    box('off');

    fname = 'CRonsetLatencyTsc1.eps';
    saveFigures(gca, fname);
    fname2 = 'CRonsetLatencyTsc1.png';
    saveFigures(gca, fname2);

    CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

    subplot(2,2,1)
    binEdges = 0:25:350;
    %CR onset time
    IkUtils.histogramIK(WT.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(WT.onsets(p.lastDay,:), binEdges) %IkUtils.histogramIK(WT.onsets(p.lastDay-1,:))
    xlim([0 350])
    title('WT onset latency')
    legend('day 1', sprintf('day%d', p.lastDay))

    subplot(2,2,2)
    %heatmap/density plot
    IkUtils.histogramIK(MUT.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(MUT.onsets(p.lastDay,:),binEdges) %IkUtils.histogramIK(MUTL7.onsets(p.lastDay-1,:)) % right now day 9 because better results...
    xlim([0 350])
    title('MUT onset latency')

    subplot(2,2,3)
    binEdges = 0:25:350;
    %CR onset time
    IkUtils.histogramIK(WTL7.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(WTL7.onsets(p.lastDay,:), binEdges) %IkUtils.histogramIK(WTL7.onsets(p.lastDay-1,:))
    xlim([0 350])
    title('WTL7 onset latency')
    legend('day 1', sprintf('day%d', p.lastDay))

    subplot(2,2,4)
    %heatmap/density plot
    IkUtils.histogramIK(MUTL7.onsets(1,:), binEdges)
    hold on
    IkUtils.histogramIK(MUTL7.onsets(p.lastDay,:),binEdges) %IkUtils.histogramIK(MUTL7.onsets(p.lastDay-1,:)) % right now day 9 because better results...
    xlim([0 350])
    title('MUTL7 onset latency')


    fname = 'CRpercentageTsc1.eps';
    saveFigures(gca, fname);
    fname2 = 'CRpercentageTsc1.png';
    saveFigures(gca, fname2);

end

function codes = getCodesMouseGroup(type)
defaultMiceList = defaultMice();
codes = [defaultMiceList([defaultMice().type] == type)];
codes = [codes.code];
end

function trial_data = trial_data_recalc(trial_data)
for f = 1:length(trial_data)
    trial_data(f).meanAll = nanmean(trial_data(f).tracesnorm);
    trial_data(f).meanCS_US = nanmean(trial_data(f).tracesnorm(trial_data(f).c_usdur~=0,:));
    trial_data(f).meanCSonly = nanmean(trial_data(f).tracesnorm(trial_data(f).c_usdur==0,:));
    trial_data(f).meanUSonly = nanmean(trial_data(f).tracesnorm(trial_data(f).c_csdur==0,:));

    %% CR amp calculation
    traces_norm_new = trial_data(f).tracesnorm;
    CSselection=traces_norm_new(trial_data(f).c_csdur==270,75:90);
    trial_data(f).CRamp=nanmax(traces_norm_new(trial_data(f).c_csdur==270,41:90)');
    baseline_2(f) = nanmean(nanmean(traces_norm_new(trial_data(f).c_csdur==0,41:90)));
    trial_data(f).CRamp5 = trial_data(f).CRamp(trial_data(f).CRamp>(1.5*baseline_2(f)));

    trial_data(f).meanCR = nanmean(trial_data(f).CRamp);
    trial_data(f).stdCR = nanstd(trial_data(f).CRamp);
    trial_data(f).meanCR5 = nanmean(trial_data(f).CRamp5);
    trial_data(f).stdCR5 = nanstd(trial_data(f).CRamp5);
end
end

function traces_new = normalizeEyelidTraces(behavior_trial_data)

baseline = nanmean(nanmean(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,1:40)));
traces_new = behavior_trial_data.eyelidpos - (nanmean(behavior_trial_data.eyelidpos(:,1:40), 2) - baseline);
baseline_new = nanmean(nanmean(traces_new(behavior_trial_data.c_csdur==0,1:40)));
stdev = nanstd(behavior_trial_data.eyelidpos(:,1:40), 0, 2);
mask_stdev = stdev < 3*nanmean(stdev);
mask_lowest = min(traces_new(:,41:end), [], 2) > (baseline_new - 0.3);
mask_highest = max(traces_new(:,141:end), [],2) < max(max(traces_new(:,90:140)));
traces_new(~mask_stdev | ~mask_lowest | ~mask_highest, :) = nan;
min_val = nanmin(traces_new(:));
max_val = nanmax(traces_new(:));
traces_new = (traces_new - min_val) / (max_val - min_val);

end

function  WT = calcMice(WTmice, p)
type = ones(1,p.nTrials);
WTtraces=nan(p.lastDay,size(WTmice,2)*p.nTrials,p.nTimeSteps);
WTcramp=nan(p.lastDay,size(WTmice,2),p.nCSplusRegTrials);
WTcramp5=nan(p.lastDay,size(WTmice,2),p.nCSplusRegTrials);
WTmean=nan(p.lastDay,p.nTimeSteps);
WTonsets=nan(p.lastDay,size(WTmice,2)*p.nTrials);

for i=1:size(WTmice,2)
    mouse = WTmice(i);
    dir_content = dir(fullfile(p.pathBehaviorDataHome, mouse));
    path_list = IkUtils.getPathFromDir(dir_content);
    for f = 1 : length(path_list)
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mouse, s, folder);
        load(fullfile(folder, filename(end)));
        trial_data(f) = behavior_trial_data;
        traces_new(f).traces = normalizeEyelidTraces(behavior_trial_data);
        trial_data(f).tracesnorm = traces_new(f).traces;
    end
    trial_data_old = trial_data;

    trial_data = trial_data_recalc(trial_data);

    if length(trial_data(f).CRamp)<220
        trial_data(f).CRamp(length(trial_data(f).CRamp)+1:220) = nan;
    end



    for j=1:p.lastDay
        WTtraces(j,(i*p.nTrials-(p.nTrials-1)):i*p.nTrials,:) = trial_data(j).tracesnorm(1:p.nTrials,:);
        WTcramp(j,i,1:min(p.nCSplusRegTrials, length(trial_data(j).CRamp))) = trial_data(j).CRamp(1:min(p.nCSplusRegTrials, length(trial_data(j).CRamp)))';
        WTcramp5(j,i,1:length(trial_data(j).CRamp5)) = trial_data(j).CRamp5';
        WTCRperc(j,i) = sum(~isnan(WTcramp5(j,i,:)))/p.nCSplusRegTrials*100;

        %         WT(i).s(j).traces = trial_data(j).tracesnorm(1:p.nTrials,:);
        %         WT(i).s(j).cramp = trial_data(j).CRamp(1:min(p.nCSplusRegTrials, length(trial_data(j).CRamp)))';
        %         WT(i).s(j).cramp5 = trial_data(j).CRamp5';
        %         WT(i).s(j).CRperc = sum(~isnan(WTcramp5(j,i,:)))/p.nCSplusRegTrials*100;
    end
end
WTCRperc(WTCRperc == 0) = nan;

type(trial_data(1).c_csdur == 0) = 0;
type(trial_data(1).c_usdur == 0) = 2;
WTtype = ones(length(WTmice), p.nTrials);
for i = 1 : length(WTmice)
    WTtype(i, :) = type;
    %     WT(i).type = type;
end

for j =1:p.lastDay
    WTcrampmean(j)=nanmean(WTcramp(j,:,:), 'all');
    WTcramp5mean(j)=nanmean(WTcramp5(j,:,:), 'all');
    WTcrampstd(j) = nanstd(WTcramp(j,:,:),0, 'all');
    WTcramp5std(j) = nanstd(WTcramp5(j,:,:),0,'all');
    WTmean(j,:) = nanmean(WTtraces(j,WTtype(1:end) == 1,:),2);
end

WTCRperdstd = std(WTCRperc,1,2);
ts=trial_data(1).tm(1,:);

%CR onset
for i = 1:p.lastDay
    for j=1:(length(WTtraces))
        if ~isnan(WTtraces(i,j,1))
            baseline = nanmean(WTtraces(i,j,1:40));
            stdbase = nanstd(WTtraces(i,j,1:40));

            if find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1)
                WTonsets(i,j) = ts(find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1));
            end
        end
    end
end

WTonsets = WTonsets(:,WTtype(1:end) ==1 |WTtype(1:end) ==2);

WT.onsets = WTonsets;
WT.traces = WTtraces;
WT.cramp = WTcramp;
WT.cramp5 = WTcramp5;
WT.CRperc = WTCRperc;
WT.type = WTtype;
WT.crampmean = WTcrampmean;
WT.cramp5mean = WTcramp5mean;
WT.crampstd = WTcrampstd;
WT.cramp5std = WTcramp5std;
WT.mean = WTmean;
WT.CRperdstd = WTCRperdstd;
WT.ts = ts;
WT.onsets = WTonsets;
end

%%
%% IK 2-11-23
function batchProcessTrialsExampleIK(mInput)
close all
debug = 0;
seeFigures = 1;
plotfigures = 1;
trainingOnly = 0;

p = IkUtils.getParams();

if nargin < 1

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to process?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );
    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);

else
    if any(contains(mInput, p.mouseNames))
        mcode = p.mouseList(find(p.mouseNames == mInput));
    else
        mcode = mInput;
    end
end

try
    dir_content = dir(fullfile(p.pathBehaviorDataHome, mcode));
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    trainingPath = "/home/i.klinkhamer/Documents/Data/behaviorDataTrainingOnlyMice/";
    dir_content = dir(fullfile(trainingPath, mcode));
    path_list = IkUtils.getPathFromDir(dir_content);
    trainingOnly = true;
end

for f = 1 : length(path_list)

    folder = path_list(f);
    [folderPath, folderName, ~] = fileparts(folder);
    if ~exist(folder, 'dir')
        status = mkdir(folder);
    end
    if ~trainingOnly
        resultsFolderPattern = "behaviorDataResults";
    else
        resultsFolderPattern = "behaviorDataResultsTrainingOnlyMice";
    end
    behaviorDataResultsDir = fullfile("/home/i.klinkhamer/Documents/Data/", resultsFolderPattern);
    behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
    if ~exist(behaviorDataResultsFolder, 'dir')
        status = mkdir(behaviorDataResultsFolder);
    end
    s = find(p.s == erase(erase(folder, fileparts(folder)),"/"));
    if isempty(s)
        s = erase(erase(folder, fileparts(folder)),"/");
    end
    stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyelidpos', file_extension = '.mat');
    stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');

    if (~exist(fullfile(folder, 'trialdata.mat'), 'file') && ~exist(fullfile(folder, stamptrial), 'file')) || debug == 1
        behavior_trial_data = processTrials(folder, 'recalibrate'); % Recalibrate eyelid
    elseif seeFigures == 1
        try
            load(fullfile(folder, stamptrial));
        catch
            load(fullfile(folder, 'trialdata.mat'));
        end
    else
        continue
    end


    if ~isempty(behavior_trial_data)

        eyelidpos = behavior_trial_data.eyelidpos;
        min_val_eye = min(eyelidpos(:));
        max_val_eye = max(eyelidpos(:));
        eyelidpos = (eyelidpos - min_val_eye) / (max_val_eye - min_val_eye);
        % Define the factor by which you want to upsample
        upsample_factor = 10;
        % Create the original time points
        original_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2));
        % Create the new time points after upsampling
        new_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2)*upsample_factor);
        % Interpolate the data
        try
            eyelidpos_upsampled = interp1(original_time, eyelidpos.', new_time, 'spline');

            min_val_eye = min(eyelidpos_upsampled(:));
            max_val_eye = max(eyelidpos_upsampled(:));
            eyelidpos = (eyelidpos_upsampled - min_val_eye) / (max_val_eye - min_val_eye);
        catch
        end

        save(fullfile(folder, stampeye), 'eyelidpos');
        if (~exist(fullfile(folder, 'trialdata.mat'), 'file') || ~exist(fullfile(folder, stamptrial), 'file')) || debug == 1
            save(fullfile(folder, stamptrial), 'behavior_trial_data');
        end

        %                 try
        first_part_folder = fileparts(fileparts(fileparts(folder)));
        middle_folder = resultsFolderPattern; %'behaviorDataResults';
        end_folder = erase(folder, fileparts(fileparts(folder)));
        new_folder = fullfile(first_part_folder,middle_folder,end_folder);
        save(fullfile(new_folder, stampeye), 'eyelidpos');

        if ~exist(fullfile(new_folder, stamptrial), 'file') || debug == 1
            save(fullfile(new_folder, stamptrial), 'behavior_trial_data');
        end
        %                 catch
        %                     disp("not saved in behavior results folder")
        %                 end


        %%
        traces = behavior_trial_data.eyelidpos;
        smoothedTraces = smootheTraces(traces);

        plotFigureEyeTraces(behavior_trial_data.tm(1,:), smoothedTraces)
        plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces)
        P = IkUtils.getParams();
        baseMin = P.tracesRanges.baseline.min;
        baseMax = P.tracesRanges.baseline.max;
        USmin = P.tracesRanges.us.min;
        USmax = P.tracesRanges.us.max;
        CSmin = P.tracesRanges.cs.min;

        baseline_min = nanmin(...
            nanmean(traces(: ,baseMin:baseMax)));

        meanUR = nanmean(nanmax(traces(:,USmin:USmax)'));

        fullClosureRange = meanUR - baseline_min;

        maskEyeMoreThanHalfOpen = nanmean(smoothedTraces(:,baseMin:baseMax)') < 0.5 * fullClosureRange;

        smoothedTraces(~maskEyeMoreThanHalfOpen,:) = nan;

        traces_new = smoothedTraces;
        traces_2baseline = traces_new - ...
            (nanmean(traces_new(:,baseMin:baseMax), 2) - ...
            baseline_min);
        % baseline_new = nanmean(nanmean(traces_new(behavior_trial_data.c_csdur==0,baseMin:baseMax)));
        baseline_new_min = nanmin(...
            nanmean(traces_2baseline(behavior_trial_data.c_csdur==0,baseMin:baseMax)));
        % stdev = nanstd(traces(:,baseMin:baseMax), 0, 2);
        % stdev = nanstd(traces(:,baseMin:baseMax), 0, 2);
        % mask_stdev = stdev < 3*nanmean(stdev);
        mask_lowest = min(traces_2baseline(:,CSmin:end), [], 2) > (baseline_new_min - 0.3);
        mask_highest = max(traces_2baseline(:,USmax+1:end), [],2) < max(max(traces_2baseline(:,USmin:USmax)));
        % traces_2baseline(~mask_lowest | ~mask_highest, :) = nan;
        traces_new(~mask_lowest | ~mask_highest, :) = nan;
        min_val = nanmin(traces_new(:));
        max_val = nanmax(traces_new(:));
        traces_norm = (traces_new - min_val) / (max_val - min_val);

        plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces_norm)

        traces_2baseline(~mask_lowest | ~mask_highest, :) = nan;
        min_val = nanmin(traces_2baseline(:));
        max_val = nanmax(traces_2baseline(:));
        traces_2baseline_norm = (traces_2baseline - min_val) / (max_val - min_val);
        plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces_2baseline_norm)
        %%

        baseline = mean(mean(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,1:40)));

        fullclosure = mean(max(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,40:95),[],2));
        traces = (behavior_trial_data.eyelidpos - baseline) ./ (fullclosure - baseline);

        traces_new = behavior_trial_data.eyelidpos - (mean(behavior_trial_data.eyelidpos(:,1:40), 2) - baseline);
        baseline_new = min(mean(traces_new(behavior_trial_data.c_csdur==0,1:40)));
        stdev = std(behavior_trial_data.eyelidpos(:,1:40), 0, 2);
        mask_stdev = stdev < 3*mean(stdev);
        mask_lowest = min(traces_new(:,41:end), [], 2) > (baseline_new - 0.3);
        mask_highest = max(traces_new(:,141:end), [],2) < max(max(traces_new(:,90:140)));
        traces_new = traces_new(mask_stdev & mask_lowest & mask_highest, :);
        min_val = min(traces_new(:));
        max_val = max(traces_new(:));
        traces_new = (traces_new - min_val) / (max_val - min_val);

        if plotfigures
            plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces_new)
        end
    end
end

end

function smoothedTraces = smootheTraces(traces)
% Parameters
windowSize = 7; % Window size for Savitzky-Golay filter
degree = 2; % Degree of polynomial for Savitzky-Golay filter

smoothedTraces = sgolayfilt(traces', degree, windowSize);

alpha = 0.05; % Significance level for Grubb's test
maxIterations = 5; % Maximum number of iterations for outlier removal
% Iterative Grubb's outlier detection
baselineSDs = std(smoothedTraces(1:40, :));
baselineMeans = mean(smoothedTraces(1:40, :));
isOutlier = false(size(baselineSDs));
for iter = 1:maxIterations
    meanSD = mean(baselineSDs);
    stdSD = std(baselineSDs);
    zScoresSD = abs((baselineSDs - meanSD) / stdSD);
    [~, idxSD] = max(zScoresSD);

    meanmean = mean(baselineMeans);
    stdMean = std(baselineMeans);
    zScoresMean = abs((baselineMeans - meanmean) / stdMean);
    [~, idxMean] = max(zScoresMean);

    zScoreThreshold = tinv(1 - alpha / (2 * numel(baselineSDs)), numel(baselineSDs) - 1);
    if zScoresSD(idxSD) > zScoreThreshold
        isOutlier(idxSD) = true;
        baselineSDs(idxSD) = [];
        baselineMeans(idxSD) = [];
    end
    if zScoresMean(idxMean) > zScoreThreshold
        isOutlier(idxMean) = true;
        baselineSDs(idxMean) = [];
        baselineMeans(idxMean) = [];
    else
        break;
    end
end
% Remove trials with unstable baselines
if any(isOutlier)
    fprintf('Removing unstable trials: %d\n', find(isOutlier));
    smoothedTraces(:, isOutlier) = nan;
else
    fprintf('No unstable trials detected.\n');
end
smoothedTraces = smoothedTraces';
end

function plotFigureEyeTraces(tm, traces)
figure; hold on
x1 = 0;
x2 = 250;
y1 = 0;
y2 = 1;
patch ...
    ( ...
    [x1 x1 x2 x2] ...
    , [y1 y2 y2 y1] ...
    , 'blue' ...
    , EdgeColor = 'none' ...
    , FaceAlpha = 0.2 ...
    )
x1 = 250;
x2 = 270;
patch ...
    ( ...
    [x1 x1 x2 x2] ...
    , [y1 y2 y2 y1] ...
    , 'green' ...
    , EdgeColor = 'none' ...
    , FaceAlpha = 0.2 ...
    )
plot(tm, traces)

title("Eyelid closure")
end

%%
eyeblinkTracesOverviewL7 = figure('Color', 'white','Position', [10 0 1000 700]);
subplot(2,2,1)
hold on
WTL7tracesPlot = WTL7.traces(:,:,:);%WTL7.traces(:,481:720,:);
plot(WT.ts(1:40), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,1:40)),'color', [0.5 0.5 0.5])
plot(WT.ts(41:90), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,41:90)),'color', [0 0 1])
plot(WT.ts(91:95), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(WT.ts(95:p.nTimeSteps), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(WT.ts, nanmean(squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: WTL7', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,2)
hold on
MUTL7tracesPlot = MUTL7.traces(:,:,:);%MUTL7.traces(:,481:720,:);
plot(WT.ts(1:40), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,1:40)),'color', [0.5 0.5 0.5])
plot(WT.ts(41:90), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,41:90)),'color', [0 0 1])
plot(WT.ts(91:95), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(WT.ts(95:p.nTimeSteps), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(WT.ts, nanmean(squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: MUTL7', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')
subplot(2,2,3)
WTL7meanWF=WTL7.mean;
WTL7meanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
WTL7meanWF(:,p.nTimeSteps)=NaN;
colormap(winter)
for i=1:p.lastDay
    color(i,1:p.nTimeSteps)=i*(length(WTL7meanWF)); %give each line a different color (instead of gradient for heigth profile)
end
h1=waterfall(WT.ts,1:p.lastDay,WTL7meanWF,color);
h1.LineWidth=2;
h1.FaceAlpha=0; %transparent faces
view([-0.5 -1.5 1]) %determines view angle
zlim=([0 1.2]);
yticks(1:p.lastDay);
xlim([-200 800]);
title(sprintf('Average traces day 1-%d: WTL7', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

subplot(2,2,4)
MUTL7meanWF=MUTL7.mean;
MUTL7meanWF(:,1)=NaN;
MUTL7meanWF(:,p.nTimeSteps)=NaN;
h2=waterfall(WT.ts,1:p.lastDay,MUTL7meanWF,color);
h2.LineWidth=2;
h2.FaceAlpha=0;
view([-0.5 -1.5 1])
zlim=([0 1.2]);
yticks(1:p.lastDay);
title(sprintf('Average traces day 1-%d: MUTL7', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

fname = sprintf('eyeblinkTracesOverview%sL7_%s.eps', Gene, IkUtils.Now);
saveFigures(eyeblinkTracesOverviewL7, fname);
fname2 = sprintf('eyeblinkTracesOverview%sL7_%s.png', Gene, IkUtils.Now);
saveFigures(eyeblinkTracesOverviewL7, fname2);eyeblinkTracesOverviewL7 = figure('Color', 'white','Position', [10 0 1000 700]);
subplot(2,2,1)
hold on
WTL7tracesPlot = WTL7.traces(:,:,:);%WTL7.traces(:,481:720,:);
plot(WT.ts(1:40), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,1:40)),'color', [0.5 0.5 0.5])
plot(WT.ts(41:90), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,41:90)),'color', [0 0 1])
plot(WT.ts(91:95), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(WT.ts(95:p.nTimeSteps), squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(WT.ts, nanmean(squeeze(WTL7tracesPlot(p.lastDay,WTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: WTL7', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,2)
hold on
MUTL7tracesPlot = MUTL7.traces(:,:,:);%MUTL7.traces(:,481:720,:);
plot(WT.ts(1:40), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,1:40)),'color', [0.5 0.5 0.5])
plot(WT.ts(41:90), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,41:90)),'color', [0 0 1])
plot(WT.ts(91:95), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(WT.ts(95:p.nTimeSteps), squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(WT.ts, nanmean(squeeze(MUTL7tracesPlot(p.lastDay,MUTL7.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: MUTL7', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')
subplot(2,2,3)
WTL7meanWF=WTL7.mean;
WTL7meanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
WTL7meanWF(:,p.nTimeSteps)=NaN;
colormap(winter)
for i=1:p.lastDay
    color(i,1:p.nTimeSteps)=i*(length(WTL7meanWF)); %give each line a different color (instead of gradient for heigth profile)
end
h1=waterfall(WT.ts,1:p.lastDay,WTL7meanWF,color);
h1.LineWidth=2;
h1.FaceAlpha=0; %transparent faces
view([-0.5 -1.5 1]) %determines view angle
zlim=([0 1.2]);
yticks(1:p.lastDay);
xlim([-200 800]);
title(sprintf('Average traces day 1-%d: WTL7', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

subplot(2,2,4)
MUTL7meanWF=MUTL7.mean;
MUTL7meanWF(:,1)=NaN;
MUTL7meanWF(:,p.nTimeSteps)=NaN;
h2=waterfall(WT.ts,1:p.lastDay,MUTL7meanWF,color);
h2.LineWidth=2;
h2.FaceAlpha=0;
view([-0.5 -1.5 1])
zlim=([0 1.2]);
yticks(1:p.lastDay);
title(sprintf('Average traces day 1-%d: MUTL7', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

fname = sprintf('eyeblinkTracesOverview%sL7_%s.eps', Gene, IkUtils.Now);
saveFigures(eyeblinkTracesOverviewL7, fname);
fname2 = sprintf('eyeblinkTracesOverview%sL7_%s.png', Gene, IkUtils.Now);
saveFigures(eyeblinkTracesOverviewL7, fname2);

%%

eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);
reg_trials_mask = repmat(WT.type(1,1:p.nTrials)==1, 1, length(WTmice));
subplot(2,2,1)
hold on
plot(ts(1:41), squeeze(WT.traces(1,WT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
plot(ts(41:91), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
plot(ts(91:95), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(ts(95:p.nTimeSteps), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(ts, nanmean(squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: WT', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,2)
hold on
MUTtracesPlot = MUT.traces(:,:,:);%MUT.traces(:,481:720,:);
plot(ts(1:41), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
plot(ts(41:91), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
plot(ts(91:95), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
plot(ts(95:p.nTimeSteps), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(ts, nanmean(squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
title(sprintf('All traces day %d: MUT', p.lastDay))
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,3)
WTmeanWF=WT.mean;
WTmeanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
WTmeanWF(:,p.nTimeSteps)=NaN;
colormap(winter)
for i=1:p.lastDay
    color(i,1:p.nTimeSteps)=i*(length(WTmeanWF)); %give each line a different color (instead of gradient for heigth profile)
end
h1=waterfall(ts,1:p.lastDay,WTmeanWF,color)
h1.LineWidth=2;
h1.FaceAlpha=0; %transparent faces
view([-0.5 -1.5 1]) %determines view angle
zlim=([0 1.2]);
yticks(1:p.lastDay);
xlim([-200 800]);
title(sprintf('Average traces day 1-%d: WT', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

subplot(2,2,4)
MUTmeanWF=MUT.mean;
MUTmeanWF(:,1)=NaN;
MUTmeanWF(:,p.nTimeSteps)=NaN;
h2=waterfall(ts,1:p.lastDay,MUTmeanWF,color)
h2.LineWidth=2;
h2.FaceAlpha=0;
view([-0.5 -1.5 1])
zlim=([0 1.2]);
yticks(1:p.lastDay);
title(sprintf('Average traces day 1-%d: MUT', p.lastDay))
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

%%

% p = IkUtils.getParams();
%     eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);
% %     reg_trials_mask = repmat(WT.type(1,1:p.nTrials)==1, 1, length(WTmice));
%     subplot(2,2,1)
%     hold on
%     plot(WT.ts(1:41), squeeze(WT.traces(1,WT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
%     plot(WT.ts(41:91), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
%     plot(WT.ts(91:95), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
%     plot(WT.ts(95:p.nTimeSteps), squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
%     plot(WT.ts, nanmean(squeeze(WT.traces(p.lastDay,WT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
%     title(sprintf('All traces day %d: %sWT', p.lastDay, localization))
%     xlabel('Time from CS (ms)')
%     ylabel('fraction eyelid closure')
%
%     subplot(2,2,2)
%     hold on
%     MUTtracesPlot = MUT.traces(:,:,:);%MUTtraces(:,481:720,:);
%     plot(WT.ts(1:41), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
%     plot(WT.ts(41:91), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,41:91)),'color', [0 0 1])
%     plot(WT.ts(91:95), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,91:95)),'color', [0 1 0])
%     plot(WT.ts(95:p.nTimeSteps), squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,95:p.nTimeSteps)),'color', [0.5 0.5 0.5])
%     plot(WT.ts, nanmean(squeeze(MUTtracesPlot(p.lastDay,MUT.type(:,1:p.nTrials)==1,:))), 'k', 'LineWidth',2);
%     title(sprintf('All traces day %d: %sMUT', p.lastDay, localization))
%     xlabel('Time from CS (ms)')
%     ylabel('fraction eyelid closure')
%
%
%     subplot(2,2,3)
%     WTmeanWF=WT.mean;
%     WTmeanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
%     WTmeanWF(:,p.nTimeSteps)=NaN;
%     colormap(winter)
%     for i=1:p.lastDay
%         color(i,1:p.nTimeSteps)=i*(length(WTmeanWF)); %give each line a different color (instead of gradient for heigth profile)
%     end
%     h1=waterfall(WT.ts,1:p.lastDay,WTmeanWF,color);
%     h1.LineWidth=2;
%     h1.FaceAlpha=0; %transparent faces
%     view([-0.5 -1.5 1]) %determines view angle
%     zlim=([0 1.2]);
%     yticks(1:p.lastDay);
%     xlim([-200 800]);
%     title(sprintf('Average traces day 1-%d: %sWT', p.lastDay, localization))
%     xlabel('Time from CS (ms)')
%     zlabel('fraction eyelid closure')
%     ylabel('training day')
%
%     subplot(2,2,4)
%     MUTmeanWF=MUT.mean;
%     MUTmeanWF(:,1)=NaN;
%     MUTmeanWF(:,p.nTimeSteps)=NaN;
%     h2=waterfall(WT.ts,1:p.lastDay,MUTmeanWF,color);
%     h2.LineWidth=2;
%     h2.FaceAlpha=0;
%     view([-0.5 -1.5 1])
%     zlim=([0 1.2]);
%     yticks(1:p.lastDay);
%     title(sprintf('Average traces day 1-%d: %sMUT', p.lastDay, localization))
%     xlabel('Time from CS (ms)')
%     zlabel('fraction eyelid closure')
%     ylabel('training day')

p = IkUtils.getParams();
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

subplot(2,2,1)
hold on
errorbar(WT.crampmean,WT.crampstd, 'b');
errorbar(MUT.crampmean,MUT.crampstd, 'r');
xticks(1:p.lastDay);
ylim([0 100]);
a = get(gca,'XTickLabel');
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR amplitude')
title('CS-US & CS-only trials')
box('off');
legend('WT', 'MUT')

subplot(2,2,2)
hold on
reshaped_CRperc_mean = reshape(nanmean(WT.CRperc,2), [size(nanmean(WT.CRperc,2), 1), size(nanmean(WT.CRperc,2), 3)]);
reshaped_CRperc_sd = reshape(nanstd(WT.CRperc,1,2), [size(nanstd(WT.CRperc,1,2), 1), size(nanstd(WT.CRperc,1,2), 3)]);

reshaped_CRperc_mean_MUT = reshape(nanmean(MUT.CRperc,2), [size(nanmean(MUT.CRperc,2), 1), size(nanmean(MUT.CRperc,2), 3)]);
reshaped_CRperc_sd_MUT = reshape(nanstd(MUT.CRperc,1,2), [size(nanstd(MUT.CRperc,1,2), 1), size(nanstd(MUT.CRperc,1,2), 3)]);
%     plot(reshaped_CRperc)
errorbar(reshaped_CRperc_mean,reshaped_CRperc_sd, 'b');
%     errorbar(nanmean(MUT.CRperc,2),nanstd(MUT.CRperc,1,2), 'r');
errorbar(reshaped_CRperc_mean_MUT ,reshaped_CRperc_sd_MUT, 'r');
xticks(1:p.lastDay);
ylim([0 100]);
a = get(gca,'XTickLabel');
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR percentage')
title('CS-US & CS-only trials')
box('off');

subplot(2,2,3)
binEdges = 0:25:350;
%CR onset time
IkUtils.histogramIK(WT.onsets(1,:), binEdges)
hold on
IkUtils.histogramIK(WT.onsets(p.lastDay,:), binEdges) %IkUtils.histogramIK(WT.onsets(p.lastDay-1,:))
xlim([0 350])
title('WT onset latency')
legend('day 1', sprintf('day%d', p.lastDay))

subplot(2,2,4)
%heatmap/density plot
IkUtils.histogramIK(MUT.onsets(1,:), binEdges)
hold on
IkUtils.histogramIK(MUT.onsets(p.lastDay,:),binEdges) %IkUtils.histogramIK(MUT.onsets(p.lastDay-1,:)) % right now day 9 because better results...
xlim([0 350])
title('MUT onset latency')

fname = sprintf('CRoverview%s%s_%s.eps', Gene, localization, IkUtils.Now);
saveFigures(gca, fname);
fname2 = sprintf('CRoverview%s%s_%s.png', Gene, localization, IkUtils.Now);
saveFigures(gca, fname2);

%%
if find(WTtraces(i,j,:)>(baseline+p.thresCRperc*stdbase) , 1)
    WTonsets(i,j) = ts(find(WTtraces(i,j,:)>(baseline+p.thresCRperc*stdbase) , 1));
end
%%
errorbar(1:numel(WT.CRperc),nanmean(WT.CRperc,2),nanstd(WT.CRperc,1,2), 'b');
errorbar(nanmean(MUT.CRperc,2),nanstd(MUT.CRperc,1,2), 'r');
%%
CR onset
for day = 1:p.lastDay
    for mouse=1:size(WTmice,2)
        for t = 1:p.nCSplusRegTrials
            WTtracesWithCS = WTtraces(:,WTtype(1:end) ==1 |WTtype(1:end) ==2, :);
            trial = p.nCSplusRegTrials * (mouse-1) + t;
            if ~isnan(WTtracesWithCS(day,trial,1)) && ~isnan(WTcramp5(day,mouse,t))
                baseline = nanmean(WTtracesWithCS(day,trial,baseMin:baseMax));
                CRPeak = nanmax(WTtracesWithCS(day,trial,CSmin:CSmax));
                threshold = baseline + (CRPeak - baseline) * p.thresCRperc / 100;
                CROnsetIndex = find(WTtracesWithCS(day,trial,CSmin:CSmax) > threshold, 1);
                if ~isempty(CROnsetIndex)
                    WTonsets(day,trial) = ts(CROnsetIndex+baseMax);
                else
                    WTonsets(day,trial) = NaN; % No CR onset found
                end
            end
        end
    end
end
%%
   
fname = sprintf('CRoverview%s%s_%s.eps', Gene, localization, IkUtils.Now);
saveas(CROverview, fullfile("/home/i.klinkhamer/Documents/Figures/", fname))
fname2 = sprintf('CRoverview%s%s_%s.png', Gene, localization, IkUtils.Now);
saveas(CROverview, fullfile("/home/i.klinkhamer/Documents/Figures/", fname2))
%%
fileExtension = '*.continuous';
if length(dir_contents) > 1
    channelPart_ = cellfun(@str2double, regexp(parts(:,2), '\d+', 'match'));
    formatPart_ = regexp(parts(:,2), '[^\d]+', 'match');
    formatPart = string(cell2mat(formatPart_{1}));

    formatPart_ = char({parts{:,:,2}});
    channelPart = str2double(string(formatPart_(:,end)));
    channelPart = channelPart(~isnan(channelPart));
    channelPart = channelPart_(channelPart_ == ch);
    if parts{:,1,3} == "2"; formatExtra = '_2'; else; formatExtra = ""; end
if parts(1,3) == "2"; formatExtra = '_2'; else; formatExtra = ""; end
else
    formatPart_ = regexp(parts(:,2), '[^\d]+', 'match');
    formatPart = string(cell2mat(formatPart_{1}));
    formatPart_ = char(parts{2});
    numbers2 = (regexp(test(:,2), '\d+', 'match'));
    numericArray = cellfun(@str2double, numbers2(:));
channelPart_ = cellfun(@str2double, regexp(parts(:,2), '\d+', 'match'));
channelPart = channelPart_(channelPart_ == ch);
    channelPart = str2double(string(formatPart_(end)));
    channelPart = channelPart(~isnan(channelPart));
    if parts{3} == "2"; formatExtra = '_2'; else; formatExtra = ""; end
if parts(1,3) == "2"; formatExtra = '_2'; else; formatExtra = ""; end
end
formatPart = string(formatPart_(1, 1:end-1));


filename = fullfile(path, sprintf("100_%s%d%s%s", formatPart, channelPart, formatExtra, fileExtension));
dir_file = dir(filename);

filepath = fullfile(dir_file.folder, dir_file.name);
%%
dataFileName = fullfile(mname, P.s(session), d);
if isfile(fullfile(which(dataFileName)))
    dataPath = fileparts(which(dataFileName));
    break
end
%%
% Source and destination directory paths
% mcode = "MI23.05990.05";
onlyTraining = ["MI23.00245.05"];
mouseListCopying = onlyTraining;
for m = 1 : length(mouseListCopying)
    mcode = mouseListCopying(m);
    fprintf("Starting copying for mouse: %s\n", mcode)
    mcodeParts = split(mcode, ".");
    if length(mcodeParts) > 1
        part1_ = char(mcodeParts(1));
        part1 = string(part1_(1)) + string(part1_(2));
        NASmcode = part1 + mcodeParts(2) + mcodeParts(3);
    else
        mcodeParts = split(mcode, "-");
        NASmcode = mcodeParts(2) + mcodeParts(3);
    end

    defaultMiceList = defaultMice();
    try
        if length(defaultMiceList([defaultMiceList.code] == mcode).dates) > 2 && (defaultMiceList([defaultMiceList.code] == mcode).dates(end) == defaultMiceList([defaultMiceList.code] == mcode).dates(end-1))
            secondEphysDay2Done = true;
        else
            secondEphysDay2Done = false;
        end
    catch
        secondEphysDay2Done = false;
    end

    sourceDir = sprintf("/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
    destDir = sprintf('/home/i.klinkhamer/Documents/Data/behaviorDataTrainingOnlyMice/%s/', mcode);

   
    secondSessionOfTheDay = 0;
   

    % List all files in the source directory
    folderList = dir(fullfile(sourceDir)); 
    folders = "";
    idx = 1;
    for i = 1 : length(folderList)
        if folderList(i).isdir && ~strcmp(folderList(i).name, '.') && ~strcmp(folderList(i).name, '..')
            folders(idx) = string(fullfile(folderList(i).folder, folderList(i).name));
            idx = idx + 1;
        end
    end

    for j = 1:length(folders)
        secondSessionOfTheDay = false;
        fileList = dir(folders(j));
        folderExtensions = "";
        for jj = 1 : length(fileList)
            [~, ~, folderExt] = fileparts(fileList(jj).name);
            folderExtensions(jj) = string(folderExt);
        end
        fileList = fileList(~[fileList.isdir] & folderExtensions == ".mat"); % You can specify the file extension you want to copy
        
        skipfile = zeros(1,length(fileList));
        lastSession = 1;
        fileSessions = [];
        for ii = 1:numel(fileList)
            parts = split(fileList(ii).name, '_');
            pattern = 's\d+\d+';
            mask = regexp(parts, pattern);
            matching_idx = find(cellfun(@(x) ~isempty(x), mask));
            if isempty(matching_idx)
                skipfile(ii) = true;
                continue
            end
            parts2 = split(parts(matching_idx), "");
            fileSessions(ii) = str2double(string(parts2{3}) + string(parts2{4}));
            lastSession = max(lastSession, str2double(string(parts2{3}) + string(parts2{4})));
        end
        fileSessions = unique(fileSessions);
        if j == 12 && secondEphysDay2Done
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession - 1);
                pattern2New = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession - 1);
                pattern2New = sprintf("s%d", lastSession);
            end
            secondSessionOfTheDay = true;
        else
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession);
                pattern2New = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession);
                pattern2New = sprintf("s%d", lastSession);
            end
        end
        try
            partsFile = split([fileList.name], '_');
        catch
            disp("You need to connect the NAS")
        end
        pattern = 's\d+\d+';
        maskPattern = regexp(partsFile, pattern);
        matching_indices = find(cellfun(@(x) ~isempty(x), maskPattern));
        newMaskPatternFiles = partsFile(matching_indices);
        newMaskFilestest = [];
        for p = 1:length(fileSessions)
            pSes = fileSessions(p);
            patternNewtest = sprintf("s0%d", pSes);
            newMaskFilestest(p).files = strcmp(newMaskPatternFiles, patternNewtest);
        end
        completed_sessions_mask = arrayfun(@(x) sum(x.files) >= 482 & sum(x.files) <= 483 | sum(x.files) >= 723, newMaskFilestest);
        if sum(completed_sessions_mask) > 1
            secondSessionOfTheDay = true;
            disp("test")
        end
        completedSessions = fileSessions(completed_sessions_mask == 1);
        newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(1)));
        newMaskFiles2 = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(end)));
        if ~(sum(newMaskFiles) >= 482 && sum(newMaskFiles) <= 482 || sum(newMaskFiles) >= 723)
            if lastSession < 10
                patternNew = sprintf("s0%d", lastSession);
            else
                patternNew = sprintf("s%d", lastSession);
            end
            newMaskFiles = ~strcmp(newMaskPatternFiles, patternNew) & ~newMaskFiles;
        end
        if length(newMaskFiles) < length(fileList)
            nExtraIdcs = length(fileList) - length(newMaskFiles);
            newMaskFiles(length(newMaskFiles)+1:length(newMaskFiles)+nExtraIdcs) = 0;
            newMaskFiles2(length(newMaskFiles2)+1:length(newMaskFiles2)+nExtraIdcs) = 0;
        end
        % Loop through each file and copy it to the destination directory
        for i = 1:numel(fileList)
   
            if newMaskFiles(i)
                [~, folderName, ~] = fileparts(folders(j));
                destFolder = fullfile(destDir, folderName);
                if ~exist(destFolder, 'dir')
                    status = mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    status = mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
            if secondSessionOfTheDay && newMaskFiles2(i)
                [~, folderName, ~] = fileparts(folders(j));
                folderNameNew = sprintf("%s2", folderName);
                destFolder = fullfile(destDir, folderNameNew);
                if ~exist(destFolder, 'dir')
                    status = mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderNameNew);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    status = mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
        end
        fprintf("folder %d done \n", j)
    end

     batchProcessTrialsExampleIK(mcode)
     status = rmdir(destDir, 's');     
end

ephysSources = [...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-02_14-02-50_1/Record Node 101/"...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-03_14-40-18_1/Record Node 101/"...
    ];
ephysDest = "/home/i.klinkhamer/Documents/Data/ephysData/";
for s = 1 : length(ephysSources)
    ephysSource = ephysSources(s);
    fileListEphys = dir(ephysSource);
    fileListEphys = fileListEphys(~[fileListEphys.isdir]);
    erasePart = fileparts(fileparts(fileparts(ephysSource)));
    ephysDestNew = fullfile(ephysDest, erase(ephysSource, erasePart));

    if ~exist(ephysDestNew, 'dir')
        status = mkdir(ephysDestNew);
    end

    for f = 1:length(fileListEphys)
        EphysSourceFile = fullfile(ephysSource, fileListEphys(f).name);
        EphysDestFile = fullfile(ephysDestNew, fileListEphys(f).name);
        if ~exist(EphysDestFile, 'file')
            copyfile(EphysSourceFile, EphysDestFile);
        end
    end
end



%%
addpath(genpath(p.pathSpikeSortingDATA));
addpath(genpath(p.pathSpikeSortingHome));
addpath(genpath(p.pathEphysDataDATA));
addpath(genpath(p.pathEphysDataHome));
addpath(genpath(p.pathBehaviorDataDATA));
addpath(genpath(p.pathBehaviorDataHome));
addpath(p.prbPath);
addpath(fullfile(path_current, "+IkUtils"));
addpath(fullfile(path_current, "DCN suite revised"));
addpath(fullfile(path_current, "Code Nynke"));

%%
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
CSmin = P.tracesRanges.cs.min;

% Calculate baselines as the mean of traces within the specified range
baseline_min = min(nanmean(traces(:, baseMin:baseMax), 2));

meanUR = nanmean(nanmax(traces(:,USmin:USmax)'));

fullClosureRange = meanUR - baseline_min;

maskEyeMoreThanHalfOpen = nanmean(smoothedTraces(:,baseMin:baseMax)') < 0.5 * fullClosureRange;

smoothedTraces(~maskEyeMoreThanHalfOpen,:) = nan;

traces_new = smoothedTraces;
traces_2baseline = traces_new - ...
    (nanmean(traces_new(:,baseMin:baseMax), 2) - ...
    baseline_min);
% baseline_new = nanmean(nanmean(traces_new(behavior_trial_data.c_csdur==0,baseMin:baseMax)));
baseline_new_min = nanmin(...
    nanmean(traces_2baseline(behavior_trial_data.c_csdur==0,baseMin:baseMax)'));
% stdev = nanstd(traces(:,baseMin:baseMax), 0, 2);
stdev = nanstd(traces_2baseline(:,baseMin:baseMax), 0, 2);
mask_stdev = stdev < 3*nanmean(stdev);
mask_lowest = min(traces_2baseline(:,CSmin:end), [], 2) > (baseline_new_min - 0.1);
mask_highest = max(traces_2baseline(:,USmax+1:end), [],2) < max(max(traces_2baseline(:,USmin:USmax)'));
% traces_2baseline(~mask_lowest | ~mask_highest, :) = nan;
traces_new(~mask_lowest | ~mask_highest |~mask_stdev, :) = nan;
min_val = nanmin(traces_new(:)');
max_val = nanmax(traces_new(:)');
traces_norm = (traces_new - min_val) / (max_val - min_val);

%                 plotFigureEyeTraces(behavior_trial_data.tm(1,:), traces_norm)

traces_2baseline(~mask_lowest | ~mask_highest|~mask_stdev, :) = nan;
min_val = nanmin(nanmean(traces_2baseline(:,baseMin:baseMax)'));
max_val = nanmax(nanmax(traces_2baseline(:,USmin:USmax)'));
traces_2baseline_norm = (traces_2baseline - min_val) / (max_val - min_val);
%%

function traces_norm = normalizeEyelidTraces(behavior_trial_data)
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
CSmin = P.tracesRanges.cs.min;

traces = behavior_trial_data.eyelidpos;

baseline_min = nanmin(...
    nanmean(traces(: ,baseMin:baseMax)));

meanUR = nanmean(nanmax(traces(:,USmin:USmax)'));

fullClosureRange = meanUR - baseline_min;

smoothedTraces = smootheTraces(traces);

maskEyeMoreThanHalfOpen = nanmean(smoothedTraces(:,baseMin:baseMax)') < 0.5 * fullClosureRange;

smoothedTraces(~maskEyeMoreThanHalfOpen,:) = nan;

traces_2baseline = smoothedTraces - ...
    (nanmean(smoothedTraces(:,baseMin:baseMax), 2) - ...
    baseline_min);
baseline_new_min = nanmin(...
    nanmean(traces_2baseline(behavior_trial_data.c_csdur==0,baseMin:baseMax)));
stdev = nanstd(traces_2baseline(:,baseMin:baseMax), 0, 2);
mask_stdev = stdev < 3*nanmean(stdev);
mask_lowest = min(traces_2baseline(:,CSmin:end), [], 2) > (baseline_new_min - 0.3);
mask_highest = max(traces_2baseline(:,USmax+1:end), [],2) < max(max(traces_2baseline(:,USmin:USmax)));
traces_2baseline(~mask_lowest | ~mask_highest | ~mask_stdev, :) = nan;
min_val = nanmin(traces_2baseline(:));
max_val = nanmax(traces_2baseline(:));
traces_norm = (traces_2baseline - min_val) / (max_val - min_val);
end
%%
 for f = 1 : length(path_list)

        folder = path_list(f);
        folder_contents = dir(fullfile(folder, 'trialdata*.mat'));
               
        s = find(P.s == erase(erase(folder, fileparts(folder)),"/"));
        if isempty(s)
            s = erase(erase(folder, fileparts(folder)),"/");
        end
       
        stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');

        if (~exist(fullfile(folder, 'trialdata.mat'), 'file') && ~exist(fullfile(folder, stamptrial), 'file')) || debug == 1
            behavior_trial_data = processTrials(folder, 'recalibrate'); % Recalibrate eyelid
        elseif seeFigures == 1
            
            load(folder_contents(1).name);            
        else
            continue
        end
%%
savefigures = true;
allMice = 0;
startMouse = 5;

if nargin < 1

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to curate?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMiceTraining() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMiceTraining() ...
        );

    mouseNames = mouseCodes + "     " + mouseNames;
    mouseNames(end+1) = "All mice";
    [~, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mouseCodes = IkUtils.promptMouseNames();
    try
        mouseCodes = mouseCodes(mname_idx);
    catch
        allMice = 1;
    end
else
    mouseCodes = IkUtils.promptMouseNames(mInput);
    mouseCodes = mInput;
end
if allMice
    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMiceTraining() ...
        );
    mouseCodes(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT") = [];
end
%%
if nargin < 1

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to process?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMiceTraining() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMiceTraining() ...
        );
    mouseNames(end+1) = "All mice";
    [~, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    try
        mouseCodes = mouseCodes(mname_idx);
    catch
        allMice = 1;
    end

else
    if any(contains(mInput, P.mouseNames))
        mouseCodes = P.mouseList(find(P.mouseNames == mInput));
    else
        mouseCodes = mInput;
    end
end
if allMice
    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMiceTraining() ...
        );
    mouseCodes(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT") = [];
end
%%
%% IK 5-4-24
function copyFromNAS(mInput)
startMouse = 1;

if nargin < 1
    mouseCodes = IkUtils.promptMouseNames();
else
    mouseCodes = IkUtils.promptMouseNames(mInput);
end
mouseCodes(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT") = [];

if length(mouseCodes) == 1
    startMouse = 1;
end

mouseListCopying = mouseCodes;
for m = startMouse : length(mouseListCopying)
    mcode = mouseListCopying(m);
    fprintf("Starting copying for mouse: %s\n", mcode)
    mcodeParts = split(mcode, ".");
    if length(mcodeParts) > 1
        part1_ = char(mcodeParts(1));
        part1 = string(part1_(1)) + string(part1_(2));
        NASmcode = part1 + mcodeParts(2) + mcodeParts(3);
    else
        mcodeParts = split(mcode, "-");
        NASmcode = mcodeParts(2) + mcodeParts(3);
    end

    defaultMiceList = defaultMiceTraining();
    try
        if length(defaultMiceList([defaultMiceList.code] == mcode).dates) > 2 && (defaultMiceList([defaultMiceList.code] == mcode).dates(end) == defaultMiceList([defaultMiceList.code] == mcode).dates(end-1))
            secondEphysDay2Done = true;
        else
            secondEphysDay2Done = false;
        end
    catch
        secondEphysDay2Done = false;
    end

    sourceDir = sprintf("/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
    destDir = sprintf('/home/i.klinkhamer/Documents/Data/behaviorDataTrainingOnlyMice/%s/', mcode);



    % List all files in the source directory
    folderList = dir(fullfile(sourceDir));
    folders = "";
    idx = 1;
    for i = 1 : length(folderList)
        if folderList(i).isdir && ~strcmp(folderList(i).name, '.') && ~strcmp(folderList(i).name, '..')
            folders(idx) = string(fullfile(folderList(i).folder, folderList(i).name));
            idx = idx + 1;
        end
    end

    for j = 1:length(folders)
        secondSessionOfTheDay = false;
        fileList = dir(folders(j));
        folderExtensions = "";
        for jj = 1 : length(fileList)
            [~, ~, folderExt] = fileparts(fileList(jj).name);
            folderExtensions(jj) = string(folderExt);
        end
        fileList = fileList(~[fileList.isdir] & folderExtensions == ".mat"); % You can specify the file extension you want to copy

        skipfile = zeros(1,length(fileList));
        lastSession = 1;
        fileSessions = [];
        for ii = 1:numel(fileList)
            parts = split(fileList(ii).name, '_');
            pattern = 's\d+\d+';
            mask = regexp(parts, pattern);
            matching_idx = find(cellfun(@(x) ~isempty(x), mask));
            if isempty(matching_idx)
                skipfile(ii) = true;
                continue
            end
            parts2 = split(parts(matching_idx), "");
            fileSessions(ii) = str2double(string(parts2{3}) + string(parts2{4}));
            lastSession = max(lastSession, str2double(string(parts2{3}) + string(parts2{4})));
        end
        fileSessions = unique(fileSessions);

        try
            partsFile = split([fileList.name], '_');
        catch
            disp("You need to connect the NAS")
        end
        pattern = 's\d+\d+';
        maskPattern = regexp(partsFile, pattern);
        matching_indices = find(cellfun(@(x) ~isempty(x), maskPattern));
        newMaskPatternFiles = partsFile(matching_indices);
        newMaskFilestest = [];
        for p = 1:length(fileSessions)
            pSes = fileSessions(p);
            if pSes < 10
                patternNewtest = sprintf("s0%d", pSes);
            else
                patternNewtest = sprintf("s%d", pSes);
            end
            newMaskFilestest(p).files = strcmp(newMaskPatternFiles, patternNewtest);
        end
        num_files_session = arrayfun(@(x) sum(x.files), newMaskFilestest);
        completed_sessions_mask = num_files_session >= 482 & num_files_session <= 483 | num_files_session >= 723;
        %         completed_sessions_mask = arrayfun(@(x) sum(x.files) >= 482 & sum(x.files) <= 483 | sum(x.files) >= 723, newMaskFilestest);
        if sum(completed_sessions_mask) > 1
            secondSessionOfTheDay = true;
            %             disp("test")
        end
        completedSessions = fileSessions(completed_sessions_mask);
        if ~isempty(completedSessions)
            if completedSessions(1) < 10
                newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(1)));
            else
                newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s%d", completedSessions(1)));
            end
            if completedSessions(end) < 10
                newMaskFiles2 = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(end)));
            else
                newMaskFiles2 = strcmp(newMaskPatternFiles, sprintf("s%d", completedSessions(end)));
            end
        else
            [~, idx] = max(num_files_session); % IK TODO: now it just takes the biggest session, but maybe add a way to combine multiple incompleted sessions.
            newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s%d", fileSessions(idx)));
        end

        if length(newMaskFiles) < length(fileList)
            nExtraIdcs = length(fileList) - length(newMaskFiles);
            newMaskFiles(length(newMaskFiles)+1:length(newMaskFiles)+nExtraIdcs) = 0;
            newMaskFiles2(length(newMaskFiles2)+1:length(newMaskFiles2)+nExtraIdcs) = 0;
        end
        % Loop through each file and copy it to the destination directory
        for i = 1:numel(fileList)

            if newMaskFiles(i)
                [~, folderName, ~] = fileparts(folders(j));
                destFolder = fullfile(destDir, folderName);
                if ~exist(destFolder, 'dir')
                    mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
            if secondSessionOfTheDay && newMaskFiles2(i)
                [~, folderName, ~] = fileparts(folders(j));
                folderNameNew = sprintf("%s2", folderName);
                destFolder = fullfile(destDir, folderNameNew);
                if ~exist(destFolder, 'dir')
                    mkdir(destFolder);
                end

                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
                behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderNameNew);
                if ~exist(behaviorDataResultsFolder, 'dir')
                    mkdir(behaviorDataResultsFolder);
                end
                sourceFolder = fullfile(sourceDir, folderName);
                sourceFile = fullfile(sourceFolder, fileList(i).name);
                destFile = fullfile(destFolder, fileList(i).name);
                if ~exist(destFile, "file")
                    copyfile(sourceFile, destFile);
                end
            end
        end
        fprintf("folder %d done \n", j)
    end

    batchProcessTrialsExampleIK(mcode)
    rmdir(destDir, 's');
end

ephysSources = [...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-02_14-02-50_1/Record Node 101/"...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-03_14-40-18_1/Record Node 101/"...
    ];
ephysDest = "/home/i.klinkhamer/Documents/Data/ephysData/";
for s = 1 : length(ephysSources)
    ephysSource = ephysSources(s);
    fileListEphys = dir(ephysSource);
    fileListEphys = fileListEphys(~[fileListEphys.isdir]);
    erasePart = fileparts(fileparts(fileparts(ephysSource)));
    ephysDestNew = fullfile(ephysDest, erase(ephysSource, erasePart));

    if ~exist(ephysDestNew, 'dir')
        mkdir(ephysDestNew);
    end

    for f = 1:length(fileListEphys)
        EphysSourceFile = fullfile(ephysSource, fileListEphys(f).name);
        EphysDestFile = fullfile(ephysDestNew, fileListEphys(f).name);
        if ~exist(EphysDestFile, 'file')
            copyfile(EphysSourceFile, EphysDestFile);
        end
    end
end


end
%%
ephysSources = [...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-02_14-02-50_1/Record Node 101/"...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-03_14-40-18_1/Record Node 101/"...
    ];
ephysDest = "/home/i.klinkhamer/Documents/Data/ephysData/";
for s = 1 : length(ephysSources)
    ephysSource = ephysSources(s);
    fileListEphys = dir(ephysSource);
    fileListEphys = fileListEphys(~[fileListEphys.isdir]);
    erasePart = fileparts(fileparts(fileparts(ephysSource)));
    ephysDestNew = fullfile(ephysDest, erase(ephysSource, erasePart));

    if ~exist(ephysDestNew, 'dir')
        mkdir(ephysDestNew);
    end

    for f = 1:length(fileListEphys)
        EphysSourceFile = fullfile(ephysSource, fileListEphys(f).name);
        EphysDestFile = fullfile(ephysDestNew, fileListEphys(f).name);
        if ~exist(EphysDestFile, 'file')
            copyfile(EphysSourceFile, EphysDestFile);
        end
    end
end
%%
function folder = sspkResultsFolder(prefixPattern)

try
    pathList = genpath(IkUtils.getParams().dirHome);
    current_file_path = which(mfilename);
    stop = 0;
    folder = [];
    loop = 0;
    tic
    while ~stop 
        loop = loop + 1;
        current_file_path = fileparts(current_file_path);
        directory = dir(fullfile(current_file_path, '**', prefixPattern + '*'));
        if ~isempty(directory)
            folder = directory.folder;
            stop = true;
        elseif length(current_file_path) <= length(char("/home/"))
            stop = true;
            %disp("Previous curation folder not found")
        end
        if loop >= 3
            stop = true;
            %disp("Previous curation folder not found")
        end
    end   
    if isempty(folder)
        folder = "";
        %disp("Previous curation folder not found")
    end
catch
    %disp("Previous curation folder not found")
    folder = "";
    return
end
% P = IkUtils.getParams(); % IK change
% %     folder = fullfile ...  % IK change
% %         ( Env.getBayesLabUserRoot() ...
% %         , "TraceExperiments" ...
% %         , "ComplexSpikeToolkit" ...
% %         , "CSpk_Suite_revised" ...
% %         , "curationResults" ...
% %         );
% folder = fullfile(P.dirHome); % IK change
end
%%
function  WT = calcMice(WTmice, P, Gene, loc, Type)
if exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)), 'file')
    load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)))
    return
end
type = ones(1,P.nTrials);
WTtraces=nan(P.lastDay,size(WTmice,2)*P.nTrials,P.nTimeSteps);
WTcramp=nan(P.lastDay,size(WTmice,2),P.nCSplusRegTrials);
WTcramp5=nan(P.lastDay,size(WTmice,2),P.nCSplusRegTrials);
WTmean=nan(P.lastDay,P.nTimeSteps);
WTonsets=nan(P.lastDay,size(WTmice,2)*P.nCSplusRegTrials);
WTCRperc=nan(P.lastDay,size(WTmice,2),P.nCSplusRegTrials);
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;

for i=1:size(WTmice,2)
    mouse = WTmice(i);
    path_list = getBehaviorPathList(mouse);
    for f = 1 : length(path_list)
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mouse, s, folder);
        try
            load(fullfile(folder, filename(end)));
        catch
            continue
        end
        if length(behavior_trial_data.session_of_day) < P.nTrials
            addNantoTrialData(mouse);
            load(fullfile(folder, filename(end)));
        end
        try
            behavior_trial_data = rmfield(behavior_trial_data, "c_iti");
            behavior_trial_data = rmfield(behavior_trial_data, "eyelidpos_original");
        catch
        end

        trial_data(f) = behavior_trial_data;

        traces_new(f).traces = normalizeEyelidTraces(behavior_trial_data);
        trial_data(f).tracesnorm = traces_new(f).traces;
    end
    trial_data = trial_data_recalc(trial_data);

    if length(trial_data(f).CRamp)<220
        trial_data(f).CRamp(length(trial_data(f).CRamp)+1:220) = nan;
    end

    for j=1:P.lastDay
        WTtraces(j,(i*P.nTrials-(P.nTrials-1)):i*P.nTrials,:) = trial_data(j).tracesnorm(1:P.nTrials,:);
        WTcramp(j,i,1:min(P.nCSplusRegTrials, length(trial_data(j).CRamp))) = trial_data(j).CRamp(1:min(P.nCSplusRegTrials, length(trial_data(j).CRamp)))';
        WTcramp5(j,i,1:length(trial_data(j).CRamp5)) = trial_data(j).CRamp5';
        WTCRperc(j,i) = sum(~isnan(WTcramp5(j,i,:)))/P.nCSplusRegTrials*100;
    end
end
WTCRperc(WTCRperc == 0) = nan;

type(trial_data(1).c_csdur == 0) = 0;
type(trial_data(1).c_usdur == 0) = 2;
WTtype = ones(length(WTmice), P.nTrials);
for i = 1 : length(WTmice)
    WTtype(i, :) = type;
end

for j =1:P.lastDay
    WTcrampmean(j)=nanmean(WTcramp(j,:,:), 'all');
    WTcramp5mean(j)=nanmean(WTcramp5(j,:,:), 'all');
    WTcrampstd(j) = nanstd(WTcramp(j,:,:),0, 'all');
    WTcramp5std(j) = nanstd(WTcramp5(j,:,:),0,'all');
    WTmean(j,:) = nanmean(WTtraces(j,WTtype(1:end) == 1,:),2);
end

WTCRperdstd = std(WTCRperc,1,2);
ts=trial_data(1).tm(1,:);


for day = 1:P.lastDay
    WTtracesWithCS = WTtraces(day, WTtype == 1 | WTtype == 2, :); % Extract relevant traces for the day
    WTtracesWithCS = reshape(WTtracesWithCS, [size(WTtracesWithCS, 2), size(WTtracesWithCS, 3)]);
    % Compute baseline and CR peak for all trials
    baseline = nanmean(WTtracesWithCS(:,baseMin:baseMax), 2);
    CRPeak = nanmax(WTtracesWithCS(:,CSmin:CSmax)');
    % Calculate threshold based on the amplitude from baseline to CR peak
    threshold = baseline + (CRPeak' - baseline) * P.thresCRonset / 100;
    
   % Find the first time point of continuous positive eyelid velocity exceeding the threshold for each trial
    for i = 1:size(WTtracesWithCS, 1)
        % Find the indices where eyelid velocity exceeds threshold
        idx = find(WTtracesWithCS(i, CSmin:CSmax) > threshold(i), 1);
        [~, CRPeakIdx] = max(WTtracesWithCS(i, CSmin:CSmax));
        diffValues = diff(WTtracesWithCS(i, CSmin+idx-1:CSmin+CRPeakIdx-1));
        % Check if the values are increasing after the onset index
        if ~isempty(idx) %&& idx < numel(WTtracesWithCS(i, CSmin:CSmax))
            if all(diffValues > 0)
                CROnsetIndices2(i) = idx;
            else
                % Select the next index as the onset index
                nextIdx = max(find(diffValues <= 0)) + idx; %#ok<*MXFND> 
                CROnsetIndices2(i) = nextIdx;
            end
        else
            CROnsetIndices2(i) = nan;
        end
        if CROnsetIndices2(i) < 6
           CROnsetIndices2(i) = nan;
        end
    end
    
    % Update WTonsets matrix
    WTonsets(day, ~isnan(CROnsetIndices2)) = ts(CROnsetIndices2(~isnan(CROnsetIndices2)) + baseMax);

end

WTonsets(isnan(WTcramp5(1,:))) = nan;

WT.onsets = WTonsets;
WT.traces = WTtraces;
WT.cramp = WTcramp;
WT.cramp5 = WTcramp5;
WT.CRperc = WTCRperc;
WT.type = WTtype;
WT.crampmean = WTcrampmean;
WT.cramp5mean = WTcramp5mean;
WT.crampstd = WTcrampstd;
WT.cramp5std = WTcramp5std;
WT.mean = WTmean;
WT.CRperdstd = WTCRperdstd;
WT.ts = ts;
WT.onsets = WTonsets;

save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("overviewTracesData_%s%s%s.mat", Gene, loc, Type)), "WT");
end

function traces_norm = normalizeEyelidTraces(behavior_trial_data)
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
CSmin = P.tracesRanges.cs.min;

traces = behavior_trial_data.eyelidpos;

baseline_min = nanmin(...
    nanmean(traces(: ,baseMin:baseMax)));

meanUR = nanmean(nanmax(traces(:,USmin:USmax)'));

fullClosureRange = meanUR - baseline_min;

smoothedTraces = smootheTraces(traces);

maskEyeMoreThanHalfOpen = nanmean(smoothedTraces(:,baseMin:baseMax)') < 0.5 * fullClosureRange;

smoothedTraces(~maskEyeMoreThanHalfOpen,:) = nan;

traces_2baseline = smoothedTraces - ...
    (nanmean(smoothedTraces(:,baseMin:baseMax), 2) - ...
    baseline_min);
baseline_new_min = nanmin(...
    nanmean(traces_2baseline(behavior_trial_data.c_csdur==0,baseMin:baseMax)));
stdev = nanstd(traces_2baseline(:,baseMin:baseMax), 0, 2);
mask_stdev = stdev < 3*nanmean(stdev);
mask_lowest = min(traces_2baseline(:,CSmin:end), [], 2) > (baseline_new_min - 0.3);
mask_highest = max(traces_2baseline(:,USmax+1:end), [],2) < max(max(traces_2baseline(:,USmin:USmax)));
traces_2baseline(~mask_lowest | ~mask_highest | ~mask_stdev, :) = nan;
min_val = nanmin(traces_2baseline(:));
max_val = nanmax(traces_2baseline(:));
traces_norm = (traces_2baseline - min_val) / (max_val - min_val);
end

function trial_data = trial_data_recalc(trial_data)
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
for f = 1:length(trial_data)
    traces_norm = trial_data(f).tracesnorm;
    trial_data(f).meanAll = nanmean(traces_norm); %#ok<*NANMEAN> 
    trial_data(f).meanCS_US = nanmean(traces_norm(trial_data(f).c_usdur~=0,:));
    trial_data(f).meanCSonly = nanmean(traces_norm(trial_data(f).c_usdur==0,:));
    trial_data(f).meanUSonly = nanmean(traces_norm(trial_data(f).c_csdur==0,:));

    %% CR amp calculation
    baseline_min(f) = nanmin(nanmean(traces_norm(:,baseMin:baseMax))); %#ok<*AGROW,*NANMIN> 
    meanUR(f) = nanmean(nanmax(traces_norm(trial_data(f).c_usdur~=0,USmin:USmax)')); %#ok<*UDIM> 
    fullBlinkRange = meanUR(f) - baseline_min(f);

    trial_data(f).CRamp= nanmax(traces_norm(trial_data(f).c_csdur==P.CSdur,CSmin:CSmax)'); %#ok<*NANMAX> 
    trial_data(f).CRamp = (trial_data(f).CRamp-baseline_min(f)) / fullBlinkRange * 100;
    trial_data(f).CRamp5 = trial_data(f).CRamp(trial_data(f).CRamp > P.thresCRperc);

    trial_data(f).meanCR = nanmean(trial_data(f).CRamp);
    trial_data(f).stdCR = nanstd(trial_data(f).CRamp); %#ok<*NANSTD> 
    trial_data(f).meanCR5 = nanmean(trial_data(f).CRamp5);
    trial_data(f).stdCR5 = nanstd(trial_data(f).CRamp5);
end
end


function smoothedTraces = smootheTraces(traces)
% Parameters
windowSize = 7; % Window size for Savitzky-Golay filter
degree = 2; % Degree of polynomial for Savitzky-Golay filter

smoothedTraces = sgolayfilt(traces', degree, windowSize);

alpha = 0.05; % Significance level for Grubb's test
maxIterations = 5; % Maximum number of iterations for outlier removal
% Iterative Grubb's outlier detection
baselineSDs = std(smoothedTraces(1:40, :));
baselineMeans = mean(smoothedTraces(1:40, :));
isOutlier = false(size(baselineSDs));
for iter = 1:maxIterations
    meanSD = mean(baselineSDs);
    stdSD = std(baselineSDs);
    zScoresSD = abs((baselineSDs - meanSD) / stdSD);
    [~, idxSD] = max(zScoresSD);

    meanmean = mean(baselineMeans);
    stdMean = std(baselineMeans);
    zScoresMean = abs((baselineMeans - meanmean) / stdMean);
    [~, idxMean] = max(zScoresMean);

    zScoreThreshold = tinv(1 - alpha / (2 * numel(baselineSDs)), numel(baselineSDs) - 1);

    if zScoresSD(idxSD) > zScoreThreshold
        isOutlier(idxSD) = true;
        baselineSDs(idxSD) = [];
        baselineMeans(idxSD) = [];
        if idxSD < idxMean
            idxMean = idxMean - 1;
        end
    end
    if zScoresMean(idxMean) > zScoreThreshold
        isOutlier(idxMean) = true;
        baselineSDs(idxMean) = [];
        baselineMeans(idxMean) = [];
    else
        break;
    end
end
% Remove trials with unstable baselines
if any(isOutlier)
    smoothedTraces(:, isOutlier) = nan;
end
smoothedTraces = smoothedTraces';
end
%%

P = IkUtils.getParams();
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;

WT_cramp = WT.cramp(1:day,:,:);
MUT_cramp = MUT.cramp(1:day,:,:);

for s = 1:size(WT_cramp,1)
    for m = 1:size(WT_cramp,2)
        traces_session = WT.traces(s,(m-1)*P.nTrials + 1 : m*P.nTrials, :);
        baseline_min_CS = nanmin(nanmean(traces_session(:,WT.type(m,:) ==0,CSmin:CSmax),3));
        if isnan(baseline_min_CS)
            baseline_min_CS = nanmean(s, baseMin:baseMax);
        end
        for t = 1:size(WT_cramp, 3)
            if isnan(WT_cramp(s,m,t))
                WT_cramp(s,m,t) = baseline_min_CS;
            end
        end
    end
end

for s = 1:size(MUT_cramp,1)
    for m = 1:size(MUT_cramp,2)
        traces_session = MUT.traces(s,(m-1)*P.nTrials + 1 : m*P.nTrials, :);
        baseline_min_CS = nanmin(nanmean(traces_session(:,MUT.type(m,:) ==0,CSmin:CSmax),3));
        if isnan(baseline_min_CS)
            baseline_min_CS = nanmean(s, baseMin:baseMax);
        end
        for t = 1:size(MUT_cramp, 3)
            if isnan(MUT_cramp(s,m,t))
                MUT_cramp(s,m,t) = baseline_min_CS;
            end
        end
    end
end


[sessionsWT, miceWT, trialsWT] = size(WT_cramp);
transposed_data_WT = permute(WT_cramp, [2, 1, 3]);
reshaped_data_WT = reshape(transposed_data_WT, [miceWT, sessionsWT*trialsWT]);
session_WT = repmat((1:sessionsWT)', size(reshaped_data_WT, 2)/sessionsWT, 1);
session_WT2 = repmat(session_WT', miceWT, 1);
groupWT = repmat({'WT'}, size(reshaped_data_WT, 2), miceWT)';
reshaped_data_WT_new = reshape(WT_cramp, [sessionsWT, miceWT*trialsWT]);
groupWT_new = repmat({'WT'}, size(reshaped_data_WT_new, 2),1);

[sessionsMUT, miceMUT, trialsMUT] = size(MUT_cramp);
transposed_data_MUT = permute(MUT_cramp, [2, 1, 3]);
reshaped_data_MUT = reshape(transposed_data_MUT, [miceMUT, sessionsMUT*trialsMUT]);
reshaped_data_MUT_new = reshape(MUT_cramp, [sessionsMUT, miceMUT*trialsMUT]);
session_MUT = repmat((1:sessionsMUT)', size(reshaped_data_MUT, 2)/sessionsMUT, 1);
session_MUT2 = repmat(session_MUT', miceMUT, 1);
groupMUT = repmat({'MUT'}, size(reshaped_data_MUT, 2), miceMUT)';
groupMUT_new = repmat({'MUT'}, size(reshaped_data_MUT_new, 2),1);

concatenated_array_new = [reshaped_data_WT; reshaped_data_MUT];
concatenated_array_new_new = [reshaped_data_WT_new, reshaped_data_MUT_new];
session_new = [session_WT2; session_MUT2];
group_new = [groupWT;groupMUT];
group_new_new = [groupWT_new;groupMUT_new];
%  group_new_test = [repmat({'WT'}, miceWT,1);repmat({'MUT'}, miceMUT, 1)];
%
%  rmData = table(concatenated_array_new, group_new, session_new, 'VariableNames', {'CRamp', 'Group', 'Session'});
%  rmData_test = table(concatenated_array_new, group_new, session_new, 'VariableNames', {'CRamp', 'Group', 'Session'});
% 
%  rm = fitrm(rmData_test, 'CRamp ~ 1', 'WithinDesign', 'Session');
 
 rmData = table(group_new(:), concatenated_array_new(:), session_new(:), 'VariableNames', {'Group', 'CRamp', 'Session'});
 rmData_new = table(group_new_new(:), concatenated_array_new_new(1,:)', concatenated_array_new_new(2,:)',...
     concatenated_array_new_new(3,:)', concatenated_array_new_new(4,:)', concatenated_array_new_new(5,:)',...
     'VariableNames', {'Group', 's1', 's2', 's3', 's4', 's5' });
Time = [1 2 3 4 5];
 % Perform the repeated measures ANOVA
rm = fitrm(rmData, 'CRamp ~ 1', 'WithinDesign', 'Session');
rm_new = fitrm(rmData_new, 's1-s5 ~ Group', 'WithinDesign', Time);


% concatenated_array = cat(2, WT_cramp, MUT_cramp);
% % data = [WT.cramp; MUT.cramp];
% group = [repmat({'WT'}, numel(WT.cramp), 1); repmat({'MUT'}, numel(MUT.cramp), 1)];
% session = repmat((1:size(WT.cramp, 1))', size(WT.cramp, 2) + size(MUT.cramp, 2), 1);
% % Create a dataset (table)
% rmData = table(data, group, session, 'VariableNames', {'CRamp', 'Group', 'Session'});
% % Perform repeated measures ANOVA
% rm = fitrm(rmData, 'CRampmean ~ 1', 'WithinDesign', 'Session');
ranovaResults = ranova(rm);
    fprintf('Repeated Measures ANOVA Results CR mean amplitude %s %s:\n', Gene, loc);
    disp(ranovaResults);

%%

% for session mean
WT_cramp_mean_sess = nanmean(WT_cramp_mean_mice,2)';
MUT_cramp_mean_sess = nanmean(MUT_cramp_mean_mice,2)';
concatenated_array_mean_sess = [WT_cramp_mean_sess; MUT_cramp_mean_sess];
group_mean_sess = [{'WT'};{'MUT'}];

rmData_mean_sess = table(group_mean_sess(:), concatenated_array_mean_sess(:,1), concatenated_array_mean_sess(:,2),...
     concatenated_array_mean_sess(:,3), concatenated_array_mean_sess(:,4), concatenated_array_mean_sess(:,5),...
     'VariableNames', {'Group', 's1', 's2', 's3', 's4', 's5' });

rm_mean_sess = fitrm(rmData_mean_sess, 's1-s5 ~ Group', 'WithinDesign', Sessions);
ranovaResults_mean_sess = ranova(rm_mean_sess);
fprintf('Repeated Measures ANOVA Results mean CR amplitude per session %s %s:\n', Gene, loc);
disp(ranovaResults_mean_sess);

%%

% for session mean
WT_CRperc_mean_sess = nanmean(WT_CRperc_mean_mice,2);
MUT_CRperc_mean_sess = nanmean(MUT_CRperc_mean_mice,2);
concatenated_array_mean_sess = [WT_CRperc_mean_sess; MUT_CRperc_mean_sess];
group_mean_sess = [repmat({'WT'}, day,1);repmat({'MUT'}, day,1)];

rmData_mean_sess = table(group_mean_sess(:), concatenated_array_mean_sess(1,:), concatenated_array_mean_sess(2,:),...
     concatenated_array_mean_sess(3,:), concatenated_array_mean_sess(4,:), concatenated_array_mean_sess(5,:),...
     'VariableNames', {'Group', 's1', 's2', 's3', 's4', 's5' });

rm_mean_sess = fitrm(rmData_mean_sess, 's1-s5 ~ Group', 'WithinDesign', Sessions);
ranovaResults_mean_sess = ranova(rm_mean_sess);
fprintf('Repeated Measures ANOVA Results mean CR percentage per session %s %s:\n', Gene, loc);
disp(ranovaResults_mean_sess);

%%

function repmAnova(WT,MUT, day, Gene, loc)
P = IkUtils.getParams();
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
Sessions = table((1:day)', 'VariableNames', {'Sessions'});

WT_CRperc = WT.CRperc(1:day,:,1);
MUT_CRperc = MUT.CRperc(1:day,:,1);


WT_CRperc_mean_mice = nanmean(WT_CRperc,3);
MUT_CRperc_mean_mice = nanmean(MUT_CRperc,3);


for s = 1:day
    meanValuesSessionsWT = nanmean(WT.mean(:,baseMin:baseMax),2);
    WT_CRperc_mean_mice(s, isnan(WT_CRperc_mean_mice(s,:))) = meanValuesSessionsWT(s);
    meanValuesSessionsMUT = nanmean(MUT.mean(:,baseMin:baseMax),2);
    MUT_CRperc_mean_mice(s, isnan(MUT_CRperc_mean_mice(s,:))) = meanValuesSessionsMUT(s);
end

variable_names = [{'Group'}, strcat('s', string(1:day))];
variable_names_cell = cellfun(@(x) {x}, variable_names);

% mean per mouse
concatenated_array_mean_mice = [WT_CRperc_mean_mice, MUT_CRperc_mean_mice]';
group_mean_mice = [repmat({'WT'}, size(WT_CRperc_mean_mice, 2),1);repmat({'MUT'}, size(MUT_CRperc_mean_mice, 2),1)];

combined_cell_mean = cell(size(concatenated_array_mean_mice,1), day+1);
combined_cell_mean(:, 1:size(group_mean_mice, 2)) = group_mean_mice; % Fill the cell array with the data
combined_cell_mean(:, size(group_mean_mice, 2)+1:size(group_mean_mice, 2)+size(concatenated_array_mean_mice, 2)) = num2cell(concatenated_array_mean_mice);

rmData_mean_mice = cell2table(combined_cell_mean,'VariableNames', variable_names_cell');
rm_mean_mice = fitrm(rmData_mean_mice, sprintf('s1-s%d ~ Group',day), 'WithinDesign', Sessions);
ranovaResults_mean_mice = ranova(rm_mean_mice);
fprintf('Repeated Measures ANOVA Results day %d mean CR amplitude per mouse %s %s:\n', day, Gene, loc);
disp(ranovaResults_mean_mice);


for s = 1:size(WT_CRperc,1)
    for m = 1:size(WT_CRperc,2)
        for t = 1:size(WT_CRperc, 3)
            if isnan(WT_CRperc(s,m,t))
                WT_CRperc(s,m,t) = 0;
            end
        end
    end
end

for s = 1:size(MUT_CRperc,1)
    for m = 1:size(MUT_CRperc,2)
       for t = 1:size(MUT_CRperc, 3)
            if isnan(MUT_CRperc(s,m,t))
                MUT_CRperc(s,m,t) = 0;
            end
        end
    end
end

[sessionsWT, miceWT, trialsWT] = size(WT_CRperc);
reshaped_data_WT_new = reshape(WT_CRperc, [sessionsWT, miceWT*trialsWT]);
groupWT = repmat({'WT'}, size(reshaped_data_WT_new, 2),1);

[sessionsMUT, miceMUT, trialsMUT] = size(MUT_CRperc);
reshaped_data_MUT = reshape(MUT_CRperc, [sessionsMUT, miceMUT*trialsMUT]);
groupMUT = repmat({'MUT'}, size(reshaped_data_MUT, 2),1);

concatenated_array = [reshaped_data_WT_new, reshaped_data_MUT]';
group = [groupWT;groupMUT];

combined_cell = cell(size(concatenated_array,1), day+1);
combined_cell(:, 1:size(group, 2)) = group; % Fill the cell array with the data
combined_cell(:, size(group, 2)+1:size(group, 2)+size(concatenated_array, 2)) = num2cell(concatenated_array);

rmData = cell2table(combined_cell,'VariableNames', variable_names_cell');
% Perform the repeated measures ANOVA
rm = fitrm(rmData, sprintf('s1-s%d ~ Group', day), 'WithinDesign', Sessions);
ranovaResults = ranova(rm);
fprintf('Repeated Measures ANOVA Results day %d CR percentage %s %s:\n', day, Gene, loc);
disp(ranovaResults);
end

%%
%% IK 8-4-24
function stamp = nameDateStampFiles(kwargs)
arguments
    kwargs.mcode = '';
    kwargs.s = '';
    kwargs.file_pattern = '';
    kwargs.file_extension = '';
end
mcode = kwargs.mcode;
s = kwargs.s;
file_pattern = kwargs.file_pattern;
file_extension = kwargs.file_extension;
mice = defaultMice();

if length(split(mcode, '-')) < 3 && length(split(mcode, '.')) < 3 && ~ismember(mcode, ["Shank2KOMUT", "Shank2KOWT"])
    charmcode = char(mcode);
    midcode = string(charmcode(3:end));
    miceNames = [mice.code];
    miceSplitNames1 = split(miceNames(1:31), '-');
    miceMidNames1 = miceSplitNames1(:,:,2);
    miceMidNames1 = erase(miceMidNames1, 'MI');
    miceEndNames1 = miceSplitNames1(:,:,3);
    listMidCodes = miceMidNames1 + miceEndNames1;
    miceSplitNames2 = split(miceNames(32:end-2), '.');
    miceMidNames2 = miceSplitNames2(:,:,2);
    miceEndNames2 = miceSplitNames2(:,:,3);
    listMidCodes = [listMidCodes, miceMidNames2 + miceEndNames2];
    mcode = mice(listMidCodes == midcode).code;
end


mouse_idx = find([mice(:).code] == mcode);
datePattern = '^\d{6,7}$'; 
caught = false;

if ischar(s) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
    date = s;
elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates(s))
    date = mice(mouse_idx).ephysdates(s);
elseif isinteger(s) && s <= IkUtils.getParams().n_sessions
    date = IkUtils.getParams().s(s);
else
    date = s;
end


if isempty(file_extension) && isempty(file_pattern) && ~isempty(mcode)
    if ischar(s) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
        date = s;
        stamp = sprintf("%s_%s", mcode, date);
    elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates(s))
        date = mice(mouse_idx).ephysdates(s);
        stamp = sprintf("%s_%s", mcode, date);
        if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
            stamp = stamp + "_2";
        end
    else
        stamp = sprintf("%s", mcode);
    end
elseif isempty(file_extension) && ~isempty(file_pattern) && ~isempty(mcode)
    if ischar(s) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
        date = s;
        stamp = sprintf("%s_%s_%s", file_pattern, mcode, date);
    elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates(s))
        date = mice(mouse_idx).ephysdates(s);
        stamp = sprintf("%s_%s_%s", file_pattern, mcode, date);
        if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
            stamp = stamp + "_2";
        end
    else
        stamp = sprintf("%s_%s", file_pattern, mcode);
    end
elseif isempty(file_pattern) && ~isempty(file_extension) && ~isempty(mcode)
    if ischar(s) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
        date = s;
        stamp = sprintf("%s_%s%s", mcode, date, file_extension);
    elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates(s))
        date = mice(mouse_idx).ephysdates(s);
        stamp = sprintf("%s_%s%s", mcode, date, file_extension);
        if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
            stamp = stamp + "_2";
        end
    else
        stamp = sprintf("%s%s", mcode, file_extension);
    end
elseif ~isempty(file_pattern) && ~isempty(file_extension) && ~isempty(mcode)
    try
        if ischar(s) || regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
            date = s;
            stamp = sprintf("%s_%s_%s%s", file_pattern, mcode, date, file_extension);
        end
    catch
        caught = true;
    end
    if ~isempty(s) && ~exist('stamp', 'var') && ~isempty(mice(mouse_idx).ephysdates(s))       %caught
        date = mice(mouse_idx).ephysdates(s);
        stamp = sprintf("%s_%s_%s%s", file_pattern, mcode, date, file_extension);
        if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
            stamp = stamp + "_2";
        end
    elseif ~exist('stamp', 'var')%caught
        stamp = sprintf("%s_%s%s", file_pattern, mcode, file_extension);
    end
end
end
%%
% IK 8-4-24
function stamp = nameDateStampFiles(kwargs)
arguments
    kwargs.mcode = '';
    kwargs.s = '';
    kwargs.file_pattern = '';
    kwargs.file_extension = '';
end
mcode = kwargs.mcode;
s = kwargs.s;
file_pattern = kwargs.file_pattern;
file_extension = kwargs.file_extension;
mice = defaultMice();

if length(split(mcode, '-')) < 3 && length(split(mcode, '.')) < 3 && ~ismember(mcode, ["Shank2KOMUT", "Shank2KOWT"])
    charmcode = char(mcode);
    midcode = string(charmcode(3:end));
    miceNames = [mice.code];
    miceSplitNames1 = split(miceNames(1:31), '-');
    miceMidNames1 = miceSplitNames1(:,:,2);
    miceMidNames1 = erase(miceMidNames1, 'MI');
    miceEndNames1 = miceSplitNames1(:,:,3);
    listMidCodes = miceMidNames1 + miceEndNames1;
    miceSplitNames2 = split(miceNames(32:end-2), '.');
    miceMidNames2 = miceSplitNames2(:,:,2);
    miceEndNames2 = miceSplitNames2(:,:,3);
    listMidCodes = [listMidCodes, miceMidNames2 + miceEndNames2];
    mcode = mice(listMidCodes == midcode).code;
end


mouse_idx = find([mice(:).code] == mcode);
caught = false;

datePattern = '^\d{6,7}$'; 
if ischar(s) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
    date = s;
elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates(s))
    date = mice(mouse_idx).ephysdates(s);
    if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
        date = date + "_2";
    end
elseif isinteger(s) && s <= IkUtils.getParams().n_sessions
    date = IkUtils.getParams().s(s);
else
    date = s;
end

if ~isempty(file_pattern)
    stamp = sprintf("%s", file_pattern);
end
if ~isempty(mcode)
    if exist("stamp", 'var')
        stamp = sprintf("%s_%s", stamp, mcode);
    else
        stamp = sprintf("%s", mcode);
    end
end
if ~isempty(date)
    if exist("stamp", 'var')
        stamp = sprintf("%s_%s", stamp, string(date));
    else
        stamp = sprintf("%s", string(date));
    end
end
if ~isempty(file_extension)
    if exist("stamp", 'var')
        stamp = sprintf("%s%s", stamp, file_extension);
    else
        stamp = sprintf("%s", file_extension);
    end
end
if isempty(file_extension) && isempty(file_pattern) && ~isempty(mcode)
    if ~isempty(date)
        stamp = sprintf("%s_%s", mcode, date);
    else
        stamp = sprintf("%s", mcode);
    end
elseif isempty(file_extension) && ~isempty(file_pattern) && ~isempty(mcode)
    if ~isempty(date)
        stamp = sprintf("%s_%s_%s", file_pattern, mcode, date);
    else
        stamp = sprintf("%s_%s", file_pattern, mcode);
    end
elseif isempty(file_pattern) && ~isempty(file_extension) && ~isempty(mcode)
    if ~isempty(date)
        stamp = sprintf("%s_%s%s", mcode, date, file_extension);
    else
        stamp = sprintf("%s%s", mcode, file_extension);
    end
elseif ~isempty(file_pattern) && ~isempty(file_extension) && ~isempty(mcode)
    if ~isempty(date)
        stamp = sprintf("%s_%s_%s%s", file_pattern, mcode, date, file_extension);
    else
        stamp = sprintf("%s_%s%s", file_pattern, mcode, file_extension);
    end
end
end
%%
% I.K. 7-9-23 initiate spike sorting with ironclust
function initSpikeSorting(kwargs)
arguments
    kwargs.setDataPath = false;
    kwargs.setProbePath = false;
end

% IkUtils.addPaths(mfilename)


% prompt = "Enter mouse index: ";
% mname = input(prompt, 's');
fprintf("\n-------------------------------------")
fprintf("\nWhich mouse would you like to curate?")
fprintf("\n-------------------------------------\n")

mouseCodes = arrayfun ...          % IK change
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
mouseNames = arrayfun ...
    ( @(mouse) string(mouse.name) ...
    , defaultMice() ...
    );
%mouseNames = [mouseNames, "All of the above"];


[~, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
mcode = mouseCodes(mname_idx);

prompt2 = "Enter session number: ";
session = input(prompt2);

app = getUImain();

P = IkUtils.getParams();
dataFileTypes = ["RAWdata*.bin", "continuous*.dat"];
while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        % If the UI figure has been closed, exit the loop
        break;
    end
    appData = evalin('base', 'appData');

    setDataPath = appData.setDataPath;
    setProbePath = appData.setProbePath; 

    for d = dataFileTypes
        [files1, folders1] = recursiveFileSearch(d, mcode, session, fileparts(P.pathEphysDataDATA));
        [files2, folders2] = recursiveFileSearch(d, mcode, session);
        [files3, folders3] = recursiveFileSearch(d, mcode, session, fileparts("/mnt/Data/Ilse/Data/"));
        [files4, folders4] = recursiveFileSearch(d, mcode, session);
        files = [files1 files2 files3 files4];
        folders = [folders1 folders2 folders3 folders4];
        if ~isempty(files)
            dataPath = folders(1);
            fileName = files(1);
            %         datafilePattern_ = split(d, '*');
            %         datafilePattern = datafilePattern_(1);
            %         datafileExtension = erase(datafilePattern_(2), "*");
            break
        end
    end
    if ~exist("dataPath", "var")
        disp("data not found")
        return
    end


    % configFileTypes = ["Config_h32_oe_bin.prm", "Config_h32_oe_dat.prm", "Config_h32_oe.prm"];
    config_file_pattern = "Config_h32_oe*.prm";
    [files, folders] = recursiveFileSearch(config_file_pattern, mcode, session);
    if isempty(files)
        stamp = nameDateStampFiles(mcode = mcode, s = session);
        copyfile("/home/i.klinkhamer/Documents/Config_h32_oe.prm", fullfile(P.pathSpikeSortingHome, mcode, P.s(session), sprintf("Config_h32_oe_%s.prm", stamp)));
        kwargs.setDataPath = true;
        [files, folders] = recursiveFileSearch(config_file_pattern, mcode, session);
    end
    for c = files
        configFileName_ = fullfile(mcode, P.s(session), c);
        if isfile(fullfile(which(configFileName_)))
            configFileName = fullfile(which(configFileName_));
            break
        end
    end
    if ~exist("configFileName", "var")
        disp("Params file not found")
        return
    end

    if setDataPath
        %     stamp = nameDateStampFiles(mcode = mcode, s = session, file_pattern = datafilePattern, file_extension = datafileExtension);
        setParamDataPath(configFileName, dataPath, fileName);
    end
    if setProbePath
        probePath = which("DBC_3.1-64-H2_IK.prb");
        %     setProbePath(configFileName, P.prbPath);
        setProbePath(configFileName, probePath);
    end

    if isfile(erase(configFileName, ".prm") + "_full.prm")
        str = "irc manual " + configFileName;
    else
        str = "irc full " + configFileName;
    end

    eval(str)

    disp("End of loop. Waiting for UI input.")
end
end

function setParamDataPath(filename, dataPath, fileName)
prm_file  = fopen(string(filename),'r');
prm_file_contents = fread(prm_file, 'char*1');
prm_file_contents_new = strrep(prm_file_contents, "vcFile = ''", sprintf("vcFile = '%s'", fullfile(dataPath, fileName)));
% prm_file_contents_new = strrep(prm_file_contents, "'RAWdata.bin'", sprintf("'%s/RAWdata_%s.bin'", dataPath, stamp));
% prm_file_contents_new = strrep(prm_file_contents_new, "'continuous.dat'", sprintf("'%s/continuous_%s.dat'", dataPath, stamp));
prm_file = fopen(filename,'w');
fprintf(prm_file, '%s', prm_file_contents_new);
fclose(prm_file);
end

function setProbePath(filename, probePath)
prb_file  = fopen(filename,'r');
prm_file_contents = fread(prb_file, 'char*1');
prm_file_contents_new = strrep(prm_file_contents, "/home/mick/Desktop/Ilse/DBC_3.1-64-H2_IK.prb", sprintf('%s', probePath));
% prm_file_contents_new = strrep(prm_file_contents_new, "'continuous.dat'", fullfile(probePath, 'continuous.dat'));
prb_file = fopen(filename,'w');
fprintf(prb_file, '%s', prm_file_contents_new);
fclose(prb_file);
end
%%
% IK 7/9/23     load data and convert to RAWdata.bin file and also create timestamps.mat and stimstaps.mat files.
function convertData(kwargs)
arguments
    kwargs.previewData = false;
    kwargs.saveConvertedData = true;
    kwargs.saveStimTimestamps = true;
    kwargs.convertAllSessions = false;
    kwargs.convertAllMice = false;
    kwargs.debug = false;
end

%IkUtils.addPaths(mfilename)

P = IkUtils.getParams();


if kwargs.convertAllMice == 0
        
        fprintf("\n-------------------------------------")
        fprintf("\nWhich mouse's data would you like to convert?")
        fprintf("\n-------------------------------------\n")
        
        mouseCodes = arrayfun ...          % IK change
            ( @(mouse) string(mouse.code) ...
            , defaultMice() ...
            );
        mouseNames = arrayfun ...
            ( @(mouse) string(mouse.name) ...
            , defaultMice() ...
            );
        %mouseNames = [mouseNames, "All of the above"];

        
        [~, mcode_idx] = IkUtils.do_prompt_select_option(mouseNames);
        mcode = mouseCodes(mcode_idx);
        mousenames = mcode;
        %mcode = mouseNames(mcode_idx);
else
    mousenames = P.mouseList;
end

if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    mice = defaultMice();
    dates = mice([mice.code] == mcode).ephysdates;
    numdates = sum(~cellfun(@isempty, dates));
    sessions = 1:numdates;
end

for mcode = mousenames
    for s = sessions
        fprintf("Mouse: %s Session: %d \n", mcode, s);
        path = getPath(mcode, s);

        if kwargs.saveConvertedData == 1
            if ~isempty(dir(fullfile(path, 'RAWdata*.bin')))
                disp("There is already a RAWdata.bin file");
%             elseif ~isempty(dir(fullfile(path, 'continuous*.dat')))
%                 disp("There is already a continuous.dat file");
            else
                saveConvertedData(path);
            end
        end

        if kwargs.previewData == 1
            previewData(path);
        end

        if kwargs.saveStimTimestamps == 1
            if ~isempty(dir(fullfile(path,'stimtimes*.mat'))) && ~isempty(dir(fullfile(path,'timestamps*.mat'))) && kwargs.debug == 0
                disp("There are already timestamp files");
            else
                saveStimTimestamps(path);
            end
        end
    end
end

end

% get path
function path = getPath(mcode, s)
P = IkUtils.getParams();
stamp = nameDateStampFiles(mcode = mcode, s = s);
try
    [~, path] = recursiveFileSearch("all_channels*.events", mcode, s);
catch
    dataFileName = fullfile(mcode, P.s(s), sprintf("all_channels_%s.events", stamp));
    path = fileparts(which(dataFileName));
end
end

% get filepath
function filepath = getFilepath(path, ch)
try
    fileParts = fileparts(path);
    [~, mcode, num] = fileparts(fileParts);
    mcode = mcode + num;
    reSaveWithStamp(mcode = mcode)
catch
end
filepattern = "100_*.continuous";
dir_contents = dir(fullfile(path, filepattern));
parts_ = split({dir_contents.name}, '_');
if size(parts_,2) > 1
    Parts = reshape(parts_, [size(parts_,2), size(parts_,3)]);
    channelPart = cellfun(@str2double, regexp(Parts(:,2), '\d+', 'match'));
    formatPart_ = regexp(Parts(:,2), '[^\d]+', 'match');
else
    Parts = parts_;
    channelPart = cellfun(@str2double, regexp(Parts(2), '\d+', 'match'));
    formatPart_ = regexp(Parts(2), '[^\d]+', 'match');
end
formatPart = string(cell2mat(formatPart_{1}));
if ~isempty(formatPart) && ch > 16
    ch = ch + 32;
end
filepath = fullfile(dir_contents(channelPart == ch).folder, dir_contents(channelPart == ch).name);
end

% convert data
function saveConvertedData(path)
for i = 1:32
    filepath = getFilepath(path, i);
    if isempty(filepath)
        disp("Channel.continuous file not found. Returning to main function.")
        return
    end
    data_ = load_open_ephys_data(filepath);
    if exist('data','var')
        data(i,:) = data_(1:length(data)); %load data
    else
        data(i,:) = data_;
    end
end
[pathmcode, folderS, ~] = fileparts(path);
[~, mcode, ~] = fileparts(pathmcode);
try
    stamp = nameDateStampFiles(mcode = mcode, s = str2double(erase(folderS, 'S')));
    fileID = fopen(fullfile(path, sprintf('RAWdata_%s.bin', stamp)),'w');
catch
    fileID = fopen(fullfile(path, 'RAWdata.bin'),'w');
end
fwrite(fileID,data, 'double');
fclose(fileID);
end

% preview data
function previewData(path)
P = getParams();
for ch = 1:P.num_channels
    filepath = getFilepath(path, ch);
    if isempty(filepath)
        continue
    end
    [~, timestamps(ch,:), ~] = load_open_ephys_data(filepath); %load data
end
try
    dir_content = dir(fullfile(path, 'data*.m'));
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

% make timestamps and stimstamps files
function saveStimTimestamps(path)
filepath = getFilepath(path, 1);
if isempty(filepath)
    disp("File not found, returning to main function.")
    return
end
[~,timestamps(1,:), ~] = load_open_ephys_data(filepath); %load data

dir_events = dir(fullfile(path, 'all_channels*.events'));
[~,timestamps2, ~] = load_open_ephys_data(fullfile(path, dir_events.name));

times = timestamps2 - timestamps(1,1);

[pathmcode, folderS, ~] = fileparts(path);
[~, mcode, mcodep2] = fileparts(pathmcode);
mcode = strcat(mcode, mcodep2);
stamp = nameDateStampFiles(mcode = mcode, s = str2double(erase(folderS, 'S')));

save(fullfile(path,sprintf('stimtimes_%s.mat', stamp)),'times')

T = timestamps(:,1);
save(fullfile(path, sprintf('timestamps_%s.mat', stamp)), 'T')
end
%%
% IK 15-3-2024
function overviewTraces()
clear
close all
app = getUImain(field = "ProcessonlyEphysmiceCheckBox", value = false);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        % If the UI figure has been closed, exit the loop
        break;
    end
    appData = evalin('base', 'appData');

    showBGlines = appData.showBGlines;
    savefigures = appData.savefigures; % Set this to false if it's taking too long to plot.
    day = appData.behaviorLastDay;
    Gene = appData.gene;
    loc = appData.loc;

    if appData.enumeration == "All mice of one gene" && appData.gene == "Shank2"
        disp("Not enough Shank2 L7 mice to do this. Please pick 'All mice of both types of a gene + location' with 'Shank2' and 'KO' and try again.")
        continue
    end

    mouseCodes = appData.mousecodes;
    mouseTypes = appData.mousetypes;
    non_valid_idcs = find(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT");
    mouseCodes(non_valid_idcs) = [];
    mouseTypes(non_valid_idcs) = [];
    if isempty(mouseCodes)
        disp("No eyelid data found. Pick another mouse.")
        continue
    end

    if appData.enumeration == "All mice" || appData.enumeration == "Single mouse"
        for m = 1:length(mouseCodes)
            close all
            dataStruct = calcMouseForOverviewTraces(mouseCodes(m)); % calculations
            plotEyeTracesFiguresSingleMouse(dataStruct, mouseCodes(m), savefigures, day) % Final figures
            plotCRoverviewSingleMouse(dataStruct, mouseCodes(m), savefigures, day)
        end
    elseif appData.enumeration == "All mice of one gene"
        dataStruct = struct;
        selectedMiceTypes = sort(unique(appData.mousetypes));

        types = struct;
        types.full = ["Shank2KOWT", "Shank2KOMUT", "Shank2L7WT", "Shank2L7MUT", "Tsc1KOWT", "Tsc1KOMUT", "Tsc1L7WT", "Tsc1L7MUT"];
        types.gene = ["Shank2", "Shank2", "Shank2", "Shank2", "Tsc1", "Tsc1", "Tsc1", "Tsc1"];
        types.loc = ["KO", "KO", "L7", "L7", "KO", "KO", "L7", "L7"];
        types.type = ["WT", "MUT", "WT", "MUT", "WT", "MUT", "WT", "MUT"];

        for t = 1:length(selectedMiceTypes)
            t_parts = split(selectedMiceTypes(t), ' ');
            mouseTypeCC = strjoin(t_parts, '');
            idx = find([types.full] == mouseTypeCC);
            Type = types.type(idx);
            Gene = types.gene(idx);
            Loc = types.loc(idx);
            dataStruct(t).type = selectedMiceTypes(t);
            mouseListType = mouseCodes(mouseTypes == selectedMiceTypes(t));
            dataStruct(t).data = calcMouseForOverviewTraces(mouseListType, Gene, Loc, Type);
        end
        if Gene == "Tsc1"
            MUTmice = getCodesMouseGroup("Tsc1KOMUT");
            WTL7mice = getCodesMouseGroup("Tsc1L7WT");
            MUTL7mice = getCodesMouseGroup("Tsc1L7MUT");
            WTmice = getCodesMouseGroup("Tsc1KOWT");
        elseif Gene == "Shank2"
            MUTmice = getCodesMouseGroup("Shank2KOMUT");
            MUTmice(MUTmice == "21-MI10159-01" | MUTmice == "21-MI10159-02"| MUTmice == "Shank2KOMUT") = [];
            WTL7mice = getCodesMouseGroup("Shank2L7WT");
            MUTL7mice = getCodesMouseGroup("Shank2L7MUT");
            WTmice = getCodesMouseGroup("Shank2KOWT");
            WTmice(WTmice == "21-MI10159-06"| WTmice == "Shank2KOWT") = [];
        end

        if twoGroups
            if loc == "KO"
                WTmice = WTmice; %#ok<*ASGSL>
                MUTmice = MUTmice;
            elseif loc == "L7"
                WTmice = WTL7mice;
                MUTmice = MUTL7mice;
            end
        else

        end

        %% calculations
        if twoGroups
            WT = calcMouseForOverviewTraces(WTmice, Gene, loc, "WT");
            MUT = calcMouseForOverviewTraces(MUTmice, Gene, loc, "MUT");
        else
            WT = calcMouseForOverviewTraces(WTmice, Gene, "KO", "WT"); %#ok<*UNRCH>
            MUT = calcMouseForOverviewTraces(MUTmice, Gene, "KO", "MUT");
            WTL7 = calcMouseForOverviewTraces(WTL7mice, Gene, "L7", "WT");
            MUTL7 = calcMouseForOverviewTraces(MUTL7mice, Gene, "L7", "MUT");
        end


        %% Final figures

        if twoGroups
            plotEyeTracesFigures(WT, MUT, Gene, loc, savefigures, day)
            plotCRoverview(WT, MUT, Gene, loc, savefigures, showBGlines, day)
        else
            plotEyeTracesFigures(WT, MUT, Gene, "KO", savefigures, day)
            plotEyeTracesFigures(WTL7, MUTL7, Gene, "L7", savefigures, day)
            plotCRoverviewAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, showBGlines, day)
            plotOnsetLatencyAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, day)
        end
    end
end
end

function appData = getUIAppData()
app = uimain();
try
    mouseNames = [[defaultMice().name], "N/A"];
catch
    IkUtils.addPaths();
    mouseNames = [[defaultMice().name], "N/A"];
end
app.MouseDropDown.Items = mouseNames;
app.StartmouseDropDown.Items = mouseNames;
uiwait(app.UIFigure)
if ~isvalid(fig)
    % If the UI figure has been closed, exit the loop
    appData = [];
    return;
end
appData = evalin('base', 'appData');
end

    function codes = getCodesMouseGroup(type)
        defaultMiceList = defaultMice();
        codes = [defaultMiceList([defaultMiceList.type] == type)];
        codes = [codes.code];
end

function plotEyeTracesFiguresSingleMouse(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotEyeTracesFiguresSingleMouse')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end

eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotEyelidTraces(Mouse, mcode, eyeblinkTracesOverview, [2,1,1], Day)
plotEyeTracesOverTime(Mouse, mcode, eyeblinkTracesOverview, [2,1,2], Day)

if savefigures
    fname = sprintf('eyeblinkTracesOverview_%s.eps', mcode);
    file = fullfile(P.figPath, "Training overview eps",fname);
    saveas(eyeblinkTracesOverview,file);
    fname2 = sprintf('eyeblinkTracesOverview_%s.png', mcode);
    file = fullfile(P.figPath, "Training overview",fname2);
    saveas(eyeblinkTracesOverview,file);
end
end

function plotCRoverviewSingleMouse(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotCRoverviewSingleMouse')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure(Mouse, mcode, CROverview, [2,2,1], day = Day);
plotCRperc(Mouse, mcode, CROverview, [2,2,2], day = Day);
plotOffset(Mouse, mcode, CROverview, [2,2,3], Day);

if savefigures
    fname = sprintf('CRoverview_%s.eps', mcode);
    file = fullfile(P.figPath, "Training overview eps",fname);
    saveas(CROverview,file);
    fname2 = sprintf('CRoverview_%s.png', mcode);
    file = fullfile(P.figPath, "Training overview",fname2);
    saveas(CROverview,file);
end
end

function plotEyeTracesFigures(WT, MUT, Gene, loc, savefigures, Day) %#ok<*DEFNU> 
if nargin < nargin('plotEyeTracesFigures')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotEyelidTraces(WT, [Gene, loc, "WT"], eyeblinkTracesOverview, [2,2,1], Day)
plotEyelidTraces(MUT, [Gene, loc, "MUT"], eyeblinkTracesOverview, [2,2,2], Day)
plotEyeTracesOverTime(WT, [Gene, loc, "WT"], eyeblinkTracesOverview, [2,2,3], Day)
plotEyeTracesOverTime(MUT, [Gene, loc, "MUT"], eyeblinkTracesOverview, [2,2,4], Day)

if savefigures
    fname = sprintf('eyeblinkTracesOverview%s%s_%s.eps', Gene, loc, IkUtils.Now);
    saveFigures(eyeblinkTracesOverview, fname);
    fname2 = sprintf('eyeblinkTracesOverview%s%s_%s.png', Gene, loc, IkUtils.Now);
    saveFigures(eyeblinkTracesOverview, fname2);
end
end

function plotCRoverviewAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, showBGlines, Day)
if nargin < nargin('plotCRoverviewAllGroups')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure([WT,MUT], [Gene, "KO"], CROverview, [2,2,1], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WT, MUT], [Gene, "KO"], CROverview, [2,2,2], showBackGroundLines = showBGlines, day = Day);
plotErrorBarFigure([WTL7,MUTL7], [Gene, "L7"], CROverview, [2,2,3], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WTL7, MUTL7], [Gene, "L7"], CROverview, [2,2,4], showBackGroundLines = showBGlines, day = Day);

if savefigures
    fname = sprintf('CRpercentage%s_%s.eps', Gene, IkUtils.Now);
    saveFigures(CROverview, fname);
    fname2 = sprintf('CRpercentage%s_%s.png', Gene, IkUtils.Now);
    saveFigures(CROverview, fname2);
end
end

function plotOnsetLatencyAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, Day)
if nargin < nargin('plotOnsetLatencyAllGroups')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
LatencyOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotOffset(WT, [Gene, "KO", "WT"], LatencyOverview, [2,2,1], Day)
plotOffset(MUT, [Gene, "KO", "MUT"], LatencyOverview, [2,2,2], Day)
plotOffset(WTL7, [Gene, "L7", "WT"], LatencyOverview, [2,2,3], Day)
plotOffset(MUTL7, [Gene, "L7", "MUT"], LatencyOverview, [2,2,4], Day)

if savefigures
    fname = sprintf('CRonsetLatency%s_%s.eps', Gene, IkUtils.Now);
    saveFigures(LatencyOverview, fname);
    fname2 = sprintf('CRonsetLatency%s_%s.png', Gene, IkUtils.Now);
    saveFigures(LatencyOverview, fname2);
end
end

function plotCRoverview(WT, MUT, Gene, loc, savefigures, showBGlines, Day)
if nargin < nargin('plotCRoverview')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure([WT,MUT], [Gene, loc], CROverview, [2,2,1], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WT, MUT], [Gene, loc], CROverview, [2,2,2], showBackGroundLines = showBGlines, day = Day);
plotOffset(WT, [Gene, loc, "WT"], CROverview, [2,2,3], Day);
plotOffset(MUT, [Gene, loc, "MUT"], CROverview, [2,2,4], Day);

if savefigures
    fname = sprintf('CRoverview%s%s_%s.eps', Gene, loc, IkUtils.Now);
    saveFigures(CROverview, fname);
    fname2 = sprintf('CRoverview%s%s_%s.png', Gene, loc, IkUtils.Now);
    saveFigures(CROverview, fname2);
end
end
%%
function [files, folders] = recursiveFileSearch(filename_pattern, mcode, s, parent_directory, skip_first)
arguments
    filename_pattern = "Config_h32_oe_*.csv";
    mcode = "20-MI19442-05"
    s = ""
    parent_directory = [IkUtils.getParams().dirHomeData, IkUtils.getParams().dirDATAData];
    skip_first = false;
end
P = IkUtils.getParams();
matching_files = [];
for d = 1 : length(parent_directory)
    pDir = parent_directory(d);
    directories = dir(pDir);
    % Recursively search for files matching the pattern in all subdirectories
    done = 0;
    if mcode == "" && s == ""
        for i = 1:numel(directories)
            if directories(i).isdir && ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..')
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                directories2 = dir(folder_path);
                for ii = 1:numel(directories2)
                    if directories2(ii).isdir && ~strcmp(directories2(ii).name, '.') && ~strcmp(directories2(ii).name, '..')
                        folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name);
                        files_in_folder = dir(fullfile(folder_path, filename_pattern));
                        matching_files = [matching_files; files_in_folder];
                        directories3 = dir(folder_path);
                        for iii = 1:numel(directories3)
                            if directories3(iii).isdir && ~strcmp(directories3(iii).name, '.') && ~strcmp(directories3(iii).name, '..')
                                folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name, directories3(iii).name);
                                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                                matching_files = [matching_files; files_in_folder];
                            end
                        end
                    end
                end
            elseif ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..') && done == 0
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(parent_directory, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                done = 1;
            end
        end
    elseif strcmp(s, "")
        for i = 1:numel(directories)
            if directories(i).isdir && ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..')
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                directories2 = dir(folder_path);
                for ii = 1:numel(directories2)
                    if directories2(ii).isdir && strcmp(directories2(ii).name, mcode) && ~strcmp(directories2(ii).name, '.') && ~strcmp(directories2(ii).name, '..')
                        folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name);
                        files_in_folder = dir(fullfile(folder_path, filename_pattern));
                        matching_files = [matching_files; files_in_folder];
                        directories3 = dir(folder_path);
                        for iii = 1:numel(directories3)
                            if directories3(iii).isdir && ~strcmp(directories3(iii).name, '.') && ~strcmp(directories3(iii).name, '..')
                                folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name, directories3(iii).name);
                                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                                matching_files = [matching_files; files_in_folder];
                            end
                        end
                    end
                end
            elseif ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..') && done == 0
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(parent_directory, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                done = 1;
            end
        end
    else
        for i = 1:numel(directories)
            if directories(i).isdir && ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..')
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                directories2 = dir(folder_path);
                for ii = 1:numel(directories2)
                    if directories2(ii).isdir && strcmp(directories2(ii).name, mcode) && ~strcmp(directories2(ii).name, '.') && ~strcmp(directories2(ii).name, '..')
                        folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name);
                        files_in_folder = dir(fullfile(folder_path, filename_pattern));
                        matching_files = [matching_files; files_in_folder];
                        directories3 = dir(folder_path);
                        for iii = 1:numel(directories3)
                            if directories3(iii).isdir && strcmp(directories3(iii).name, P.s(s)) && ~strcmp(directories3(iii).name, '.') && ~strcmp(directories3(iii).name, '..')
                                folder_path = fullfile(parent_directory, directories(i).name, directories2(ii).name, directories3(iii).name);
                                files_in_folder = dir(fullfile(folder_path, filename_pattern));
                                matching_files = [matching_files; files_in_folder];
                            end
                        end
                    end
                end
            elseif ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..') && done == 0
                folder_path = fullfile(parent_directory, directories(i).name);
                files_in_folder = dir(fullfile(parent_directory, filename_pattern));
                matching_files = [matching_files; files_in_folder];
                done = 1;
            end
        end
    end


end
try
    files = string({matching_files.name});
    folders = string({matching_files.folder});
catch
    files = [];
    folders = [];
end
% filenametest = recursiveFileSearch("StructEphysData*.mat", mname, s);
end

function [matching_files, directories_new] = searchFiles(matching_files, directories)
for i = 1:numel(directories)
    if directories(i).isdir && ~strcmp(directories(i).name, '.') && ~strcmp(directories(i).name, '..')
        folder_path = fullfile(parent_directory, directories(i).name);
        files_in_folder = dir(fullfile(folder_path, filename_pattern));
        matching_files = [matching_files; files_in_folder];
        directories_new = dir(folder_path);
    end
end
end
%%
% Recursively search for files matching the pattern in all subdirectories
for i = 1:numel(directories)
    [matching_files, directories2] = searchFiles(matching_files, directories(i), filename_pattern);
    for ii = 1:numel(directories2)
        [matching_files, directories3] = searchFiles(matching_files, directories2(ii), filename_pattern);
        for iii = 1:numel(directories3)
            matching_files = searchFiles(matching_files, directories3(iii), filename_pattern);
        end
    end
end
%%
% IK 25-4-24       Recursively look for files containing a certain pattern in their filename.
function [files, folders] = recursiveFileSearch(filename_pattern, mcode, s, parent_directory, max_depth)
arguments
    filename_pattern = ''; %"Config_h32_oe_*.csv";
    mcode = ''; %"20-MI19442-05"
    s = []; %1
    parent_directory = [IkUtils.getParams().dirHomeData, IkUtils.getParams().dirDATAData];
    max_depth = 3;
end
P = IkUtils.getParams();
matching_files = [];
if isempty(filename_pattern)
    disp("Didn't enter filename pattern")
    files = [];
    folders = [];
    return
end
for d = 1 : length(parent_directory)
    pDir = parent_directory(d);
    directories = dir(pDir);
    files = directories(~[directories.isdir]);
    directories = directories([directories.isdir] & ~contains(string({directories.name}), [".", ".."]));

    matching_files = searchFilesRecursively(matching_files, directories, filename_pattern, max_depth);  % Recursively search for files matching the pattern in all subdirectories

    partsFilenamePattern = split(filename_pattern, "*");
    for f = 1 : numel(files)
        if all(arrayfun(@(p) contains(files(f).name, p), partsFilenamePattern))
            matching_files = [matching_files, files(f)];
        end
    end
    if ~isempty(mcode)
        mcodeMask = contains(string({matching_files.name}), mcode) | contains(string({matching_files.folder}), mcode);
        if any(mcodeMask)
            matching_files = matching_files(mcodeMask);
        end
    end
    if ~isempty(s)
        if isinteger(s) && s <= P.n_sessions
            sFolderMask = contains(string({matching_files.name}), P.s(s)) | contains(string({matching_files.folder}), P.s(s));
            if ~isempty(mcode)
                mice = defaultMice();
                ephysDates = mice([defaultMice().code] == mcode).ephysdates;
                date = ephysDates(s);
                sEphysMask = arrayfun(@(date) contains(string({matching_files.name}), date) | contains(string({matching_files.folder}), date), ephysDates, 'UniformOutput', false);
                sEphysMask = any(cat(3, sEphysMask{:}), 3);

                parts = split(ephysDates, "-"); ephysDatesJoined = join(parts, "", 3);
                ephysDatesJoined(ephysDatesJoined.strlength > 6) = extractAfter(ephysDatesJoined(ephysDatesJoined.strlength > 6), 2);
                sJoinedMask = arrayfun(@(date) contains(string({matching_files.name}), date) | contains(string({matching_files.folder}), date), ephysDatesJoined, 'UniformOutput', false);
                sJoinedMask = any(cat(3, sJoinedMask{:}), 3);

                fullMask = sFolderMask | sEphysMask | sJoinedMask;
            else
                fullMask = sFolderMask;
            end
        else
            s = string(s); % Convert into string in case it wasn't already a string.
            
        end
        if any(fullMask)
            matching_files = matching_files(fullMask);
        end
    end
end

if ~isempty(matching_files)
    files = string({matching_files.name});
    folders = string({matching_files.folder});
else
    files = [];
    folders = [];
end

end

function [matching_files, directories_new] = searchFiles(matching_files, directory, filename_pattern)
folder_path = fullfile(directory.folder, directory.name);
matching_files_in_folder = dir(fullfile(folder_path, filename_pattern));
matching_files = [matching_files; matching_files_in_folder];
directories_new = dir(folder_path);
directories_new = directories_new([directories_new.isdir] & ~contains(string({directories_new.name}), [".", ".."]));
end

function matching_files = searchFilesRecursively(matching_files, directories, filename_pattern, depth)
    if depth == 0
        return;
    end
    
    for i = 1:numel(directories)
        [matching_files, next_directories] = searchFiles(matching_files, directories(i), filename_pattern);
        matching_files = searchFilesRecursively(matching_files, next_directories, filename_pattern, depth - 1);
    end
end
%%
% I.K. 25/10/2023   Taking Spikesorting data outcome and putting it into a nice structure to use.
function unitConstruction(kwargs)
arguments
    kwargs.convertAllSessions = true;
    kwargs.convertAllMice = false;
    kwargs.startMouse = 1;
end

% addpath('/media/mick/DATA/Ilse/convertedData/')
%IkUtils.addPaths(mfilename)

P = IkUtils.getParams();
if kwargs.convertAllMice == 0

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to curate?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...          % IK change
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );
    %mouseNames = [mouseNames, "All of the above"];

    ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
    mouseNames = mouseNames(ephysMiceMask);
    mouseCodes = mouseCodes(ephysMiceMask);

    [~, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);
    mousecodes = mcode;
    %mcode = mouseNames(mname_idx);
else
    mousecodes = P.mouseList(kwargs.startMouse:end);
end
% if kwargs.convertAllMice == 0
%     prompt = "Enter mouse index: ";
%     mousenames = string(input(prompt, 's'));
%     if any(contains(mousenames, P.mouseNames))
%         mousenames = P.mouseList(find(P.mouseNames == mousenames));
%     end
% else
%     mousenames = P.mouseList;
%     mousenames = mousenames(kwargs.startMouse:end);
% end



if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    try
        mice = defaultMice();
        dates = mice([mice.code] == mcode).ephysdates;
        numdates = sum(~cellfun(@isempty, dates));
        sessions = 1:numdates;
    catch
        sessions = 1:P.n_sessions;
    end
end

for mcode = mousecodes
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mcode, P.mouseNames(find(P.mouseList == mcode)), s);

        fileTypes = ["Config_h32_oe_bin.csv", "Config_h32_oe_dat.csv", "Config_h32_oe.csv"];
        filePattern = "Config_h32_oe*.csv";
        [files, folders] = recursiveFileSearch(filePattern, mcode, s);
        for c = files
            FileName_ = fullfile(mcode, P.s(s), c);
            if isfile(fullfile(which(FileName_)))
                filename = fullfile(which(FileName_));
                break
            end
        end      

        if ~exist("filename", "var")
            disp(".csv file not found")
            continue
        end
        
        [filesStim, foldersStim] = recursiveFileSearch('stimtimes*.mat', mcode, s);
        if isempty(filesStim)
            [filesStim, foldersStim] = recursiveFileSearch('stimtimes*.mat', mcode, s, P.dirDATAData);
        end
%         stimFileName_ = fullfile(mcode, P.s(s), 'stimtimes.mat');
        stimFileName_ = fullfile(foldersStim, filesStim);
        if isempty(filesStim) && s <= 2
            disp("No stimtimes found, this is a Rhymdata mouse.")
            continue
        elseif ~isfile(fullfile(which(stimFileName_))) && s >= 3
            fprintf("Session %d not found \n", s)
            continue
        end
        stimFileName = fullfile(which(stimFileName_));
        try
            cell_spk = import_jrc_csv(strcat(filename));
        catch
            unit = [];
            path = fileparts(which(FileName_));
            stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'StructEphysData', file_extension = ".mat");
            save(fullfile(path,stamp),'unit')
            continue
        end
        n_neurons = length(cell_spk) - 1;
        if n_neurons == 0
            disp("There are no neurons in this session")
        end

        S = load(strcat(erase(filename, ".csv"),"_jrc.mat"));

        if S.S_clu.nClu ~= n_neurons
            disp("Number of neurons in csv file does not match jrc file")
            continue
        end

%         %% loading the ITI's to determine the US only trials
%         eyeVidFileName_ = fullfile(mcode, P.s(s), 'trialdata.mat');
%         if ~isfile(fullfile(which(eyeVidFileName_))) && s <= 2
%             disp("No eye vid trial data found.")
%         end
%         eyeVidFileName = fullfile(which(eyeVidFileName_));
% 
%         trialdata = load(eyeVidFileName);
%         try
%             ITIs = trialdata.behavior_trial_data.c_iti;
%         catch
%             disp("ITIs not found in file. Try remaking the file.")
%         end


        unit.neuron = [];
        bins = -P.t_raster_edges:0.0005:P.t_raster_edges;
        for i = 1 : n_neurons
            unit.neuron(i).spiketimes = cell_spk(i+1).t;
            unit.neuron(i).irc_trueID = cell_spk(i+1).nr;
            unit.neuron(i).channels = cell_spk(i+1).ch;
            unit.neuron(i).bins_cs = bins;
            unit.neuron(i).bins_us = bins;
            unit.irc_trueID(i) = cell_spk(i+1).nr;
        end

        unit.full_spiketimes = cell_spk(1).t;
        unit.full_clusters = S.S_clu.viClu;
        unit.bins_cs = bins;
        unit.bins_us = bins;
        unit.mcode = mcode;
        unit.name = P.mouseNames(find(P.mouseList == mcode));
        unit.session = s;

        stimtimes = load(stimFileName);
        stimtimes = stimtimes.times;
        stimtimes = [stimtimes zeros(length(stimtimes),1)];

        if length(stimtimes) <= 440
            stimtimes(1:2:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(21:22:end,3) = 1; % IK: CS only trials.
            stimtimes(22:22:end,3) = 1; % IK: CS only trials.

            stimtimesCS = stimtimes(1:2:end,1); % Take only the odd stimtimes, because those are the times that the light is turning ON.
%             if exist('ITIs','var') % Finding the stimtimes of when the CS would be during the US only trials.
%                 stimtimesUS = stimtimesCS(1:11:end) - ITIs(2:12:end);
%                 stimtimesUS_;
%             end
            stimtimesCS(:,2) = 1;
            stimtimesCS(11:11:end,2) = 0; % IK: CS only trials.
            stimtimesCS_sorted_per_condition = sortrows(stimtimesCS, 2, 'ascend');

            CS_only_trial_idcs_ON = 11:11:length(stimtimesCS);
            trialsON = 1:length(stimtimesCS);
            unit.sorted_per_condition_trials = [CS_only_trial_idcs_ON  trialsON(~ismember(trialsON,CS_only_trial_idcs_ON))];
        else
            stimtimes(3:4:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(43:44:end,3) = 1; % IK: CS only trials.
            stimtimes(44:44:end,3) = 1; % IK: CS only trials.
            stimtimes(1:44:end, 4) = 1; % IK: US only trials.
            stimtimes(2:44:end, 4) = 1; % IK: US only trials.

            stimtimesCS_ = stimtimes(3:4:end,1); % times that the light is turning ON.

            stimtimesPuffidcs_ = [1:44:length(stimtimes) 4:4:length(stimtimes)];
            wrongidcs = 44:44:length(stimtimes);
            stimtimesPuffidcs = sort(stimtimesPuffidcs_(~ismember(stimtimesPuffidcs_, wrongidcs)));
            stimtimesPuff = stimtimes(stimtimesPuffidcs, 1); % times of the puff   

            stimtimes_CS_during_US_only_trials_ = stimtimesPuff(1:11:end); % US only trials
            stimtimes_CS_during_US_only_trials = stimtimes_CS_during_US_only_trials_ - (IkUtils.getParams().isi + IkUtils.getParams().digital_US_delay)/1000;
            stimtimesCS = sortrows([stimtimesCS_; stimtimes_CS_during_US_only_trials]);
           
            
            CS_only_trial_idcs = 12:12:length(stimtimesCS);
            US_only_trial_idcs = 1:12:length(stimtimesCS);
            control_trials = [CS_only_trial_idcs US_only_trial_idcs];
            normal_trials = 1:length(stimtimesCS);
            normal_trials(control_trials) = [];
            
            unit.sorted_per_condition_trials = [CS_only_trial_idcs normal_trials US_only_trial_idcs];

            stimtimesCS(CS_only_trial_idcs, 2) = 0;
            stimtimesCS(US_only_trial_idcs, 2) = 2;
            stimtimesCS(normal_trials, 2) = 1;
            stimtimesCS_sorted_per_condition = sortrows(stimtimesCS, 2, "ascend");
            
        end

        unit.stimtimesCS_sorted_per_condition = stimtimesCS_sorted_per_condition;
        unit.stimtimesCS = stimtimesCS;

        unit.stimtimes = stimtimes;

        n_trials = length(stimtimesCS);

        % extract spiketimes and stimtimes
        for n = 1 : n_neurons
            unit.neuron(n).trial_numbers = zeros(length(unit.neuron(n).spiketimes),1);
            unit.neuron(n).trial_numbers_sorted = zeros(length(unit.neuron(n).spiketimes),1);
            for trial = 1 : n_trials
                unit.neuron(n).trial_numbers(unit.neuron(n).spiketimes > (stimtimesCS(trial) - P.t_raster_edges) &...
                    unit.neuron(n).spiketimes < (stimtimesCS(trial) + P.t_raster_edges)) =...
                    trial;
                unit.neuron(n).trial_numbers_sorted(unit.neuron(n).spiketimes > (stimtimesCS_sorted_per_condition(trial) - P.t_raster_edges) &...
                    unit.neuron(n).spiketimes < (stimtimesCS_sorted_per_condition(trial) + P.t_raster_edges)) =...
                    trial;
            end
            mask_trials = unit.neuron(n).trial_numbers > 0;
            unit.neuron(n).RasterXY_cs = zeros(2, sum(mask_trials)*3);
            unit.neuron(n).RasterXY_cs(2, 1:3:end) = unit.neuron(n).trial_numbers(mask_trials);
            unit.neuron(n).RasterXY_cs(2, 2:3:end) = unit.neuron(n).trial_numbers(mask_trials) + P.raster_spike_height;
            unit.neuron(n).RasterXY_cs(2, 3:3:end) = nan;
            unit.neuron(n).RasterXY_cs(1, 3:3:end) = nan;

            mask_trials = unit.neuron(n).trial_numbers_sorted > 0;
            unit.neuron(n).RasterXY_cs_sorted = zeros(2, sum(mask_trials)*3);
            unit.neuron(n).RasterXY_cs_sorted(2, 1:3:end) = unit.neuron(n).trial_numbers_sorted(mask_trials);
            unit.neuron(n).RasterXY_cs_sorted(2, 2:3:end) = unit.neuron(n).trial_numbers_sorted(mask_trials) + P.raster_spike_height;
            unit.neuron(n).RasterXY_cs_sorted(2, 3:3:end) = nan;
            unit.neuron(n).RasterXY_cs_sorted(1, 3:3:end) = nan;
        end

        for n = 1 : n_neurons
            for trial = 1 : n_trials
                unit.neuron(n).trial_spikes_cs(trial).trial = unit.neuron(n).spiketimes(unit.neuron(n).trial_numbers == trial) - stimtimesCS(trial);
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial) = unit.neuron(n).trial_spikes_cs(trial).trial;
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial + P.raster_spike_height) = unit.neuron(n).trial_spikes_cs(trial).trial;

                unit.neuron(n).trial_spikes_cs_sorted(trial).trial = unit.neuron(n).spiketimes(unit.neuron(n).trial_numbers_sorted == trial) - stimtimesCS_sorted_per_condition(trial);
                unit.neuron(n).RasterXY_cs_sorted(1,unit.neuron(n).RasterXY_cs_sorted(2,:) ==  trial) = unit.neuron(n).trial_spikes_cs_sorted(trial).trial;
                unit.neuron(n).RasterXY_cs_sorted(1,unit.neuron(n).RasterXY_cs_sorted(2,:) ==  trial + P.raster_spike_height) = unit.neuron(n).trial_spikes_cs_sorted(trial).trial;

                unit.neuron(n).trial_spikes_us(trial).trial = unit.neuron(n).trial_spikes_cs(trial).trial + P.t_US_offset;
                unit.neuron(n).trial_spikes_us_sorted(trial).trial = unit.neuron(n).trial_spikes_cs_sorted(trial).trial + P.t_US_offset;
            end
            unit.neuron(n).RasterXY_us(1,:) = unit.neuron(n).RasterXY_cs(1,:) + P.t_US_offset;
            unit.neuron(n).RasterXY_us_sorted(1,:) = unit.neuron(n).RasterXY_cs_sorted(1,:) + P.t_US_offset;
            unit.neuron(n).RasterXY_us(2,:) = unit.neuron(n).RasterXY_cs(2,:);
            unit.neuron(n).RasterXY_us_sorted(2,:) = unit.neuron(n).RasterXY_cs_sorted(2,:);

            unit.neuron(n).RasterXY_us_filtered = rasterFilter(unit.neuron(n).RasterXY_us);
            unit.neuron(n).RasterXY_us_sorted_filtered = rasterFilter(unit.neuron(n).RasterXY_us_sorted);

            unit.neuron(n).RasterXY_cs_filtered = rasterFilter(unit.neuron(n).RasterXY_cs);
            unit.neuron(n).RasterXY_cs_sorted_filtered = rasterFilter(unit.neuron(n).RasterXY_cs_sorted);

            stimtimesCS_sorted = unit.stimtimesCS_sorted_per_condition;
            RasterXY_spikes_cs = unit.neuron(n).RasterXY_cs_sorted(1,1:3:end);
            RasterXY_trials_cs = unit.neuron(n).RasterXY_cs_sorted(2,1:3:end);

            RasterXY_spikes_us = unit.neuron(n).RasterXY_us_sorted(1,1:3:end);
            RasterXY_trials_us = unit.neuron(n).RasterXY_us_sorted(2,1:3:end);

            for c = 1 : length(unique(stimtimesCS_sorted(:,2)))
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
                spike_counts_cs = histcounts(RasterXY_spikes_cs(ismember(RasterXY_trials_cs, trials)),unit.bins_cs);
                spike_counts_cs = spike_counts_cs/length(trials);
                spike_rates_cs = spike_counts_cs/0.0005;

                spike_counts_us = histcounts(RasterXY_spikes_us(ismember(RasterXY_trials_us, trials)),unit.bins_us);
                spike_counts_us = spike_counts_us/length(trials);
                spike_rates_us = spike_counts_us/0.0005;
                % create filter
                sigma1 = 50; % pick sigma value for the gaussian
                gaussFilter1 = gausswin(6*sigma1 + 1)';
                gaussFilter1 = gaussFilter1 / sum(gaussFilter1); % normalize
                filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter1, 'same');
                filtered_spike_rates_us = conv(spike_rates_us, gaussFilter1, 'same');
                % do an extra gauss filter to smoothen out the lines
                sigma2 = 10; % pick sigma value for the gaussian
                gaussFilter2 = gausswin(6*sigma2 + 1)';
                gaussFilter2 = gaussFilter2 / sum(gaussFilter2); % normalize
                filtered_spike_rates_cs = conv(filtered_spike_rates_cs, gaussFilter2, 'same');
                filtered_spike_rates_us = conv(filtered_spike_rates_us, gaussFilter2, 'same');

                unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates_cs;
                unit.neuron(n).psth_us_reset(c, :) = filtered_spike_rates_us;
            end


            RasterXY_spikes_cs_filtered = unit.neuron(n).RasterXY_cs_sorted_filtered(1,1:3:end);
            RasterXY_trials_cs_filtered = unit.neuron(n).RasterXY_cs_sorted_filtered(2,1:3:end);

            RasterXY_spikes_us_filtered = unit.neuron(n).RasterXY_us_sorted_filtered(1,1:3:end);
            RasterXY_trials_us_filtered = unit.neuron(n).RasterXY_us_sorted_filtered(2,1:3:end);


            for c = 1 : length(unique(stimtimesCS_sorted(:,2)))
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
                spike_counts_cs = histcounts(RasterXY_spikes_cs_filtered(ismember(RasterXY_trials_cs_filtered, trials)),unit.bins_cs);
%                 spike_counts_cs = spike_counts_cs/length(trials);
                spike_counts_cs = spike_counts_cs/sum(ismember(trials, RasterXY_trials_cs_filtered));
                spike_rates_cs = spike_counts_cs/0.0005;

                spike_counts_us = histcounts(RasterXY_spikes_us_filtered(ismember(RasterXY_trials_us_filtered, trials)),unit.bins_us);
                spike_counts_us = spike_counts_us/sum(ismember(trials, RasterXY_trials_cs_filtered));
                spike_rates_us = spike_counts_us/0.0005;
                % create filter
                sigma1 = 50; % pick sigma value for the gaussian
                gaussFilter1 = gausswin(6*sigma1 + 1)';
                gaussFilter1 = gaussFilter1 / sum(gaussFilter1); % normalize
                filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter1, 'same');
                filtered_spike_rates_us = conv(spike_rates_us, gaussFilter1, 'same');
                % do an extra gauss filter to smoothen out the lines
                sigma2 = 10; % pick sigma value for the gaussian
                gaussFilter2 = gausswin(6*sigma2 + 1)';
                gaussFilter2 = gaussFilter2 / sum(gaussFilter2); % normalize
                filtered_spike_rates_cs = conv(filtered_spike_rates_cs, gaussFilter2, 'same');
                filtered_spike_rates_us = conv(filtered_spike_rates_us, gaussFilter2, 'same');

                unit.neuron(n).psth_cs_reset_filtered(c, :) = filtered_spike_rates_cs;
                unit.neuron(n).psth_us_reset_filtered(c, :) = filtered_spike_rates_us;
            end
        end

        path = fileparts(which(FileName_));
        stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'StructEphysData', file_extension = ".mat");
        save(fullfile(path,stamp),'unit')
    end
end

end

function rasterFiltered = rasterFilter(raster, range)
if nargin < 2
    range.min = min(raster(1,:));
    range.max = max(raster(1,:));
end
binWidth = IkUtils.getParams().BinW;
nTrials = max(raster(2,:)) - IkUtils.getParams().raster_spike_height;
binEdges = range.min:binWidth:range.max;
filterMask = zeros(1,length(raster));
for t = 1 : nTrials
    spikeTimesTrial_ = raster(1,raster(2,:) == t);
    rangeMask = ...
        spikeTimesTrial_ >= range.min ...
        & spikeTimesTrial_ <= range.max;
    spikeTimesTrial = spikeTimesTrial_(rangeMask);
    amplitudes = histcounts(spikeTimesTrial, binEdges);
    baseline = mean(amplitudes);
    %rate = mean(amplitudes/nTrials/binWidth);
    baserate = baseline/binWidth;
    % nSpikesTrial(t) = length(find(raster(2,:) == t));
    if baserate > IkUtils.getParams().Sspk_base_filter %nSpikesTrial(t) > 0
        filterMask(raster(2,:) == t) = 1;
        filterMask(find(raster(2,:) == t) + 1) = 1;
        filterMask(find(raster(2,:) == t) + 2) = 1;
%     else
%         fprintf("Trial %d baserate too low: %d\n", t, baserate);
    end
end
rasterFiltered(1,:) = raster(1,find(filterMask));
rasterFiltered(2,:) = raster(2,find(filterMask));
rasterFiltered(3,:) = NaN;
rasterFiltered(3, rasterFiltered(2,:) <= 20.9) = 0;
rasterFiltered(3, rasterFiltered(2,:) > 20.9 & rasterFiltered(2,:) <= 220.9) = 1;
rasterFiltered(3, rasterFiltered(2,:) > 220.9) = 2;
end

%%
% IK 8-4-24
function params = getParams(kwargs)
arguments
    % General settings for K-means clustering and filtering
    kwargs.nClusters (1,1) double = 4; % IK note: check    % IK removed complex
    kwargs.Sspk_base_filter (1,1) double = 20;
    kwargs.Sspk_max_filter (1,1) double = 200; % IK note: this might be too low. Check
    
    % General settings for histogram plots
    kwargs.histogramModulationThreshold_sspks = 5; % 5 Hz taken from ten Brinke et al., 2017 % IK added. Check if this value is correct later.
    kwargs.minimumModulationTimeThreshold_sspks = 10; % 10 ms taken from ten Brinke et al., 2017 % IK added.
    kwargs.sspkRanges = struct ...
        ( cs = struct ...
            ( min = 0.050 ... % IK changed from 0.055
            , max = 0.250 ...
            ) ...
        , us = struct ...
            ( min = 0.250... % IK change %0.05 ...
            , max = 0.350... %0.150... % IK change. %0.070 ... % I.e. 50 milliseconds after the US event
            ) ...
        );
    
    
    % Trace event times
    kwargs.delayEventsCs = [0, 0.25];  % IK change
    kwargs.delayEventsUs = [-0.25, 0]; % IK change
    
    
    % Settings for histogram bin width
    kwargs.BinW = 0.005; % 5 ms
    kwargs.sd_thres = 3/4 * std(rand(1e6,1)); % std for uniformly distributed spiking events % IK note: check for simple
    
    % Params for cross-correlation
    kwargs.cc_range = [-0.2, 0.2];
    kwargs.cc_BinW = 0.002;
    kwargs.cc_minBin = 101;
    kwargs.cc_maxBin = 113;
    kwargs.cc_nBins = 13;  
    
    % Params for contamination
    kwargs.cThreshold = 10;
    kwargs.cBinW = 5;
    kwargs.cBinLimits = [0,1200];
    
    % Params for boxplot figures
    kwargs.bpTicks = 1:2;   % IK note: Change for DCN
    kwargs.bpTickLabels = {'Delta', 'Uniform'}; % IK note: change for DCN
    kwargs.bpFontSize = 12;
    
end

    kwargs.psthRanges = struct ...
        ( cs_full = struct ...
            ( min = -0.25 ...
            , max = 0.75... %0.75 ... % IK change
            ) ...
        , cs_zoomed = ...
            struct ...
                ( min = min(-0.02   , kwargs.sspkRanges.cs.min) ...
                , max = max(0.5 , kwargs.sspkRanges.cs.max) ... % IK change
                ) ...
        , us_full = ...
            struct ...
                ( min = -0.75 ...
                , max = 0.25 ...
                ) ...
        , us_zoomed = ...
            struct ...
                ( min = -.02 ...
                , max = 0.200 ...
                ) ...
        );
% overviewTraces
kwargs.lastDay = 10;
kwargs.nTimeSteps = 200;
kwargs.nTrials = 240;
kwargs.nCSplusRegTrials = 220;

kwargs.CSdur = 270;
kwargs.USdur = 20;

kwargs.thresCRperc = 5;
kwargs.thresCRonset = 5;

kwargs.tracesRanges = struct ...
        ( baseline = struct ...
            ( min = 1 ... 
            , max = 40 ...
            ) ...
        , cs = struct ...
            ( min = 41 ... 
            , max = 90 ...
            ) ...
        , us = struct ...
            ( min = 91 ... 
            , max = 140 ...  
            ) ...
         , after = struct ...
            ( min = 141 ...
            , max = 200 ...
            ) ...
        );

% kwargs.complex_edges = [0:kwargs.BinW:0.2]; % this no longer works since we
% have multiple complex spike search ranges (cs, prior, us).

% general
kwargs.n_sessions = 3;
kwargs.s = ["S1", "S2", "S3"];
kwargs.n_conditions = 3; % There are 2 conditions, CS-US paired and CS-only
kwargs.conditions = [0 1 2]; % For CS-only, 0, for CS-US paired 1, for US-only 2;
kwargs.condition_types = ["CS only", "CS-US combined", "US only"];
kwargs.eventTimes = [0 250];
% Paths
kwargs.dirDATA = "/mnt/Data/Ilse/";
kwargs.dirHome = "/home/i.klinkhamer/Documents/";
kwargs.dirHomeData = "/home/i.klinkhamer/Documents/Data/";
kwargs.dirDATAData = "/mnt/Data/Ilse/Data/";
kwargs.dirHomeCode = "/home/i.klinkhamer/Documents/ephys_code/";
% kwargs.pathSpikeSortingDATA = "/media/mick/DATA/Ilse/spikeSortingUnits/";
% kwargs.pathSpikeSortingHome = "/home/mick/Desktop/Ilse/spikeSortingUnits/";
% kwargs.pathEphysDataDATA = "/media/mick/DATA/Ilse/Data/convertedData/";
% kwargs.pathEphysDataHome = "/home/mick/Desktop/Ilse/convertedData/";
% kwargs.pathBehaviorDataDATA = "/media/mick/DATA/Ilse/Data/behaviorData/";
% kwargs.prbPath = '/home/mick/Desktop/Ilse/DBC_3.1-64-H2_IK.prb';
kwargs.pathSpikeSortingDATA = "/mnt/Data/Ilse/spikeSortingUnits/";
kwargs.pathSpikeSortingHome = "/home/i.klinkhamer/Documents/Data/spikeSortingUnits/";
kwargs.pathEphysDataDATA = "/mnt/Data/Ilse/Data/convertedData/";
kwargs.pathEphysDataHome = "/home/i.klinkhamer/Documents/Data/ephysData/";
kwargs.pathBehaviorDataDATA = "/mnt/Data/Ilse/Data/behaviorData/";
kwargs.pathBehaviorDataHome = "/home/i.klinkhamer/Documents/Data/behaviorData/";
kwargs.figPath = "/home/i.klinkhamer/Documents/Figures/";
kwargs.prbPath = '/home/i.klinkhamer/Documents/DBC_3.1-64-H2_IK.prb';
% Nynke's variables
kwargs.t_pre = 500;
kwargs.t_post = 250; %for CR detection
kwargs.dur = 2800;
kwargs.t_psth = 1000;
kwargs.bin_psth = 150;
kwargs.bin_corr = 250;
kwargs.isi = 250;
kwargs.digital_US_delay_ms = 12; % US appears 12 ms earlier in the digitally recorded US times than it appears in real life.
kwargs.digital_US_delay = 0.012;
kwargs.CS = 1;
kwargs.US = 3;
% psth and raster variables
kwargs.t_pre_trial = 0.500;
kwargs.t_post_trial = 2;
kwargs.t_raster_edges = 3;
kwargs.t_US_offset = 0.250;
kwargs.raster_spike_height = 0.900;
kwargs.histogramHeight = 15;
kwargs.BinW = 0.010;
kwargs.n_trials = 240;
kwargs.t_US = 0.250;
kwargs.spikeTimesBinW = 0.0005;
% convert data variables
kwargs.A = 100; %ID number of module % IK: This is the number that is in front of all .continuous files.
kwargs.num_channels = 32;
% mouselist
kwargs.mouseList = ["20-MI19154-08", "20-MI19308-01", "20-MI19601-03", '20-MI19442-03', '20-MI19442-05', '20-MI19442-06', '20-MI19442-08', '21-MI10159-01', '21-MI10159-02', '21-MI10159-06', '21-MI10532-03', '21-MI10532-07', '21-MI16091-04', '21-MI16091-03', '21-MI16183-03', '21-MI16183-05', '22-MI10447-09', '22-MI10447-05', '22-MI10008-08', '22-MI10447-06', '22-MI11756-06', '22-MI11756-07', '22-MI12417-08', '22-MI12410-07', '22-MI13134-01', '22-MI13134-03', '22-MI13989-04','22-MI13989-05', '22-MI13989-07', '22-MI14020-03', '22-MI14020-07', '22-MI14020-08', 'MI22.02438.04', 'MI22.02438.05', 'MI23.00102.01', 'MI23.00102.03', 'MI22.01127.02', 'MI22.01127.03', 'MI22.01127.04', 'MI22.01125.07', 'MI22.02178.02', 'MI22.02178.01', 'MI22.02178.04', 'MI22.02178.03', 'MI22.02170.01', 'MI23.00244.04', 'MI23.01047.03', 'MI23.01712.07', 'MI23.01753.03', 'MI23.01412.04', 'MI23.01753.02', 'MI23.02063.04', 'MI23.02063.05', 'MI23.02436.01', 'MI23.02436.03', 'MI23.02436.04', 'MI23.03123.02', 'MI23.03123.03', 'MI23.03123.04', 'MI23.03374.01', 'MI23.03374.03', 'MI23.05990.02', 'MI23.05990.03', 'MI23.05990.04', 'MI23.05990.05', 'MI23.05990.08'];
kwargs.mouseNames = ["Aurora", 'Bacchus', 'Cupid', 'Diana.1', 'Diana.2', 'Diana.3', 'Diana.4', 'Epona.1', 'Epona.2', 'Epona.3', 'Fortuna.1', 'Fortuna.2', 'Genius.2', 'Genius.1', 'Hercules.1', 'Hercules.2', 'Invidia.3', 'Invidia.1', 'Juno', 'Invidia.2', 'Luna.1', 'Luna.2', 'Mars', 'Neptune', 'Orcus.1', 'Orcus.2', 'Pluto.1', 'Pluto.2', 'Pluto.3', 'Quiritis.1', 'Quiritis.2', 'Quiritis.3', 'Roma.1', 'Roma.2', 'Saturn.1', 'Saturn.2', 'Trivia.1', 'Trivia.2', 'Trivia.3', 'Trivia.4', 'Vulcan.2', 'Vulcan.1', 'Vulcan.4', 'Vulcan.3', 'Ares', 'Bia', 'Chronos', 'Demeter', 'Eros.2', 'Fury', 'Eros.1', 'Gaia.1', 'Gaia.2', 'Hades.1', 'Hades.2', 'Hades.3', 'Icarus.1', 'Icarus.2', 'Icarus.3', 'Jason.1', 'Jason.2', 'Kratos.1', 'Kratos.2', 'Kratos.3', 'Kratos.4', 'Kratos.5'];
kwargs.Shank2MUT = ["Diana.1", 'Diana.2', 'Diana.4', 'Epona.1', 'Epona.2', 'Fortuna.2', 'Genius.2', 'Genius.1', 'Hercules.2', 'Invidia.3', 'Juno', 'Luna.1', 'Orcus.1', 'Pluto.2', 'Quiritis.2', 'Quiritis.3', 'Saturn.2', 'Bia', 'Fury'];
kwargs.Shank2WT = ["Diana.3", 'Epona.3', 'Fortuna.1', 'Hercules.1', 'Invidia.1', 'Invidia.2', 'Luna.2', 'Orcus.2', 'Pluto.1', 'Quiritis.1', 'Saturn.1', 'Chronos', 'Demeter'];

types = struct;
types.full = ["Shank2KOWT", "Shank2KOMUT", "Shank2L7WT", "Shank2L7MUT", "Tsc1KOWT", "Tsc1KOMUT", "Tsc1L7WT", "Tsc1L7MUT"];
types.gene = ["Shank2", "Shank2", "Shank2", "Shank2", "Tsc1", "Tsc1", "Tsc1", "Tsc1"];
types.loc = ["KO", "KO", "L7", "L7", "KO", "KO", "L7", "L7"];
types.type = ["WT", "MUT", "WT", "MUT", "WT", "MUT", "WT", "MUT"];
kwargs.types = types;

params = kwargs;

end
%%
% I.K. 25/10/2023   Taking Spikesorting data outcome and putting it into a nice structure to use.
function unitConstruction(kwargs)
arguments
    kwargs.convertAllSessions = true;
    kwargs.convertAllMice = false;
    kwargs.startMouse = 1;
end
P = IkUtils.getParams();
if kwargs.convertAllMice == 0

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to construct a data struct for?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = [defaultMice().code];
    mouseNames = [defaultMice().name];

    ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
    mouseNames = mouseNames(ephysMiceMask);
    mouseCodes = mouseCodes(ephysMiceMask);

%     for c = 1 : length(mouseCodes)
%         choices(c) = sprintf("%s\t\t%s",mouseNames(c), mouseCodes(c));
%     end

    choices = mouseNames;

    [~, mname_idx] = IkUtils.do_prompt_select_option(choices);
    mcode = mouseCodes(mname_idx);
    mousecodes = mcode;
else
    mousecodes = P.mouseList(kwargs.startMouse:end);
end

if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    try
        mice = defaultMice();
        dates = mice([mice.code] == mcode).ephysdates;
        numdates = sum(~cellfun(@isempty, dates));
        sessions = 1:numdates;
    catch
        sessions = 1:P.n_sessions;
    end
end

for mcode = mousecodes
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mcode, mouseNames(mouseCodes == mcode), s);

        [S, cell_spk, stimtimes] = loadData(mcode, s);
        if isempty(S) || isempty(cell_spk) || isempty(stimtimes)
            continue
        end

        unit.neuron = [];
        bins = -P.t_raster_edges : P.spikeTimesBinW : P.t_raster_edges;
        n_neurons = length(cell_spk) - 1;
        for i = 1 : n_neurons
            unit.neuron(i).spiketimes = cell_spk(i+1).t;
            unit.neuron(i).irc_trueID = cell_spk(i+1).nr;
            unit.neuron(i).channels = cell_spk(i+1).ch;
            unit.neuron(i).bins_cs = bins;
            unit.neuron(i).bins_us = bins;
            unit.irc_trueID(i) = cell_spk(i+1).nr;
        end

        unit.full_spiketimes = cell_spk(1).t;
        unit.full_clusters = S.S_clu.viClu;
        unit.bins_cs = bins;
        unit.bins_us = bins;
        unit.mcode = mcode;
        unit.name = mouseNames(mouseCodes == mcode);
        unit.session = s;

        diff_stim_times = diff(stimtimes);
        if all(diff_stim_times(2:2:end, 1) > 10) % If this is the case than only the stimtimes of the light have been recorded.
            stimtimes(1:2:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(21:22:end,3) = 1; % IK: CS only trials.
            stimtimes(22:22:end,3) = 1; % IK: CS only trials.

            stimtimesCS = stimtimes(1:2:end,1); % Take only the odd stimtimes, because those are the times that the light is turning ON.
            stimtimesCS(:,2) = 1;
            stimtimesCS(11:11:end,2) = 0; % IK: CS only trials.
            stimtimesCS_sorted_per_condition = sortrows(stimtimesCS, 2, 'ascend');

            CS_only_trial_idcs_ON = 11:11:length(stimtimesCS);
            trialsON = 1:length(stimtimesCS);
            sorted_per_condition_trials = [CS_only_trial_idcs_ON  trialsON(~ismember(trialsON,CS_only_trial_idcs_ON))];
        else % if this is the case than both the stimtimes of the the light and puff have been recorded.
            stim_times.stimtimes = stimtimes;
            stim_times.lightON = stimtimes(3:4:end);
            stim_times.lightOFF = stimtimes(4:4:end);
            stim_times.CSonlyLightON = stimtimes(43:44:end);
            stim_times.CSonlyLightOFF = stimtimes(44:44:end);

            stimtimesPuffidcs_ = [1:44:length(stimtimes) 4:4:length(stimtimes)];
            notPuffidcs = 44:44:length(stimtimes);
            stimtimesPuffidcs = sort(stimtimesPuffidcs_(~ismember(stimtimesPuffidcs_, notPuffidcs)));
            stimtimesPuff = stimtimes(stimtimesPuffidcs); % times of the puff

            stim_times.puffStart = stimtimesPuff + P.digital_US_delay;
            stim_times.puffEnd = stim_times.puffStart + (P.USdur/1000);
            stim_times.USonlyPuffStart = stimtimes(1:44:end) + P.digital_US_delay;
            stim_times.USonlyPuffEnd = stimtimes(2:44:end) + P.digital_US_delay;
            
            stimtimes(3:4:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(43:44:end,3) = 1; % IK: CS only trials.
            stimtimes(44:44:end,3) = 1; % IK: CS only trials.
            stimtimes(1:44:end, 4) = 1; % IK: US only trials.
            stimtimes(2:44:end, 4) = 1; % IK: US only trials.
            
            stimtimes_CS_during_US_only_trials = stim_times.USonlyPuffStart - (P.isi / 1000);
            stimtimesCS = sortrows([stim_times.lightON; stimtimes_CS_during_US_only_trials]);           

            stim_times.stimtimesCS = stimtimesCS;
            
            CS_only_trial_idcs = 12:12:length(stimtimesCS);
            US_only_trial_idcs = 1:12:length(stimtimesCS);
            normal_trial_idcs = setdiff(1:length(stimtimesCS), [CS_only_trial_idcs, US_only_trial_idcs]);
            
            sorted_per_condition_trials = [CS_only_trial_idcs, normal_trial_idcs, US_only_trial_idcs];

            stimtimesCS(CS_only_trial_idcs, 2) = 0;
            stimtimesCS(US_only_trial_idcs, 2) = 2;
            stimtimesCS(normal_trial_idcs, 2) = 1;
            stimtimesCS_sorted_per_condition = sortrows(stimtimesCS, 2, "ascend");

            stimtimesCS_sorted_test = stimtimesCS(sorted_per_condition_trials);
            
        end

        unit.stimtimesCS_sorted_per_condition = stimtimesCS_sorted_per_condition;
        unit.stimtimesCS = stimtimesCS;

        unit.stimtimes = stimtimes;

        n_trials = length(stimtimesCS);

        % extract spiketimes and stimtimes
        for n = 1 : n_neurons
            trial_numbers = zeros(length(unit.neuron(n).spiketimes),1);
            trial_numbers_sorted = zeros(length(unit.neuron(n).spiketimes),1);
            for trial = 1 : n_trials
                trial_numbers(unit.neuron(n).spiketimes > (stimtimesCS(trial) - P.t_raster_edges) &...
                    unit.neuron(n).spiketimes < (stimtimesCS(trial) + P.t_raster_edges)) =...
                    trial;
                trial_numbers_sorted(unit.neuron(n).spiketimes > (stimtimesCS_sorted_per_condition(trial) - P.t_raster_edges) &...
                    unit.neuron(n).spiketimes < (stimtimesCS_sorted_per_condition(trial) + P.t_raster_edges)) =...
                    trial;
            end          

            mask_trials_sorted = trial_numbers_sorted > 0;
            unit.neuron(n).RasterXY_cs_sorted = zeros(2, sum(mask_trials_sorted)*3);
            unit.neuron(n).RasterXY_cs_sorted(2, 1:3:end) = trial_numbers_sorted(mask_trials_sorted);
            unit.neuron(n).RasterXY_cs_sorted(2, 2:3:end) = trial_numbers_sorted(mask_trials_sorted) + P.raster_spike_height;
            unit.neuron(n).RasterXY_cs_sorted(1:2, 3:3:end) = nan;            

            mask_trials = trial_numbers > 0;
            unit.neuron(n).RasterXY_cs = unit.neuron(n).RasterXY_cs_sorted;
            unit.neuron(n).RasterXY_cs(2, 1:3:end) = trial_numbers(mask_trials);
            unit.neuron(n).RasterXY_cs(2, 2:3:end) = trial_numbers(mask_trials) + P.raster_spike_height;
        end

        for n = 1 : n_neurons
            for trial = 1 : n_trials
                unit.neuron(n).trial_spikes_cs(trial).trial = unit.neuron(n).spiketimes(trial_numbers == trial) - stimtimesCS(trial);
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial) = unit.neuron(n).trial_spikes_cs(trial).trial;
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial + P.raster_spike_height) = unit.neuron(n).trial_spikes_cs(trial).trial;

                unit.neuron(n).trial_spikes_cs_sorted(trial).trial = unit.neuron(n).spiketimes(trial_numbers_sorted == trial) - stimtimesCS_sorted_per_condition(trial);
                unit.neuron(n).RasterXY_cs_sorted(1,unit.neuron(n).RasterXY_cs_sorted(2,:) ==  trial) = unit.neuron(n).trial_spikes_cs_sorted(trial).trial;
                unit.neuron(n).RasterXY_cs_sorted(1,unit.neuron(n).RasterXY_cs_sorted(2,:) ==  trial + P.raster_spike_height) = unit.neuron(n).trial_spikes_cs_sorted(trial).trial;

                unit.neuron(n).trial_spikes_us(trial).trial = unit.neuron(n).trial_spikes_cs(trial).trial + P.t_US_offset;
                unit.neuron(n).trial_spikes_us_sorted(trial).trial = unit.neuron(n).trial_spikes_cs_sorted(trial).trial + P.t_US_offset;
            end
            unit.neuron(n).RasterXY_us(1,:) = unit.neuron(n).RasterXY_cs(1,:) + P.t_US_offset;
            unit.neuron(n).RasterXY_us_sorted(1,:) = unit.neuron(n).RasterXY_cs_sorted(1,:) + P.t_US_offset;
            unit.neuron(n).RasterXY_us(2,:) = unit.neuron(n).RasterXY_cs(2,:);
            unit.neuron(n).RasterXY_us_sorted(2,:) = unit.neuron(n).RasterXY_cs_sorted(2,:);

            unit.neuron(n).RasterXY_us_filtered = rasterFilter(unit.neuron(n).RasterXY_us);
            unit.neuron(n).RasterXY_us_sorted_filtered = rasterFilter(unit.neuron(n).RasterXY_us_sorted);

            unit.neuron(n).RasterXY_cs_filtered = rasterFilter(unit.neuron(n).RasterXY_cs);
            unit.neuron(n).RasterXY_cs_sorted_filtered = rasterFilter(unit.neuron(n).RasterXY_cs_sorted);

            stimtimesCS_sorted = unit.stimtimesCS_sorted_per_condition;
            RasterXY_spikes_cs = unit.neuron(n).RasterXY_cs_sorted(1,1:3:end);
            RasterXY_trials_cs = unit.neuron(n).RasterXY_cs_sorted(2,1:3:end);

            RasterXY_spikes_us = unit.neuron(n).RasterXY_us_sorted(1,1:3:end);
            RasterXY_trials_us = unit.neuron(n).RasterXY_us_sorted(2,1:3:end);

            for c = 1 : length(unique(stimtimesCS_sorted(:,2)))
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
                spike_counts_cs = histcounts(RasterXY_spikes_cs(ismember(RasterXY_trials_cs, trials)),unit.bins_cs);
                spike_counts_cs = spike_counts_cs/length(trials);
                spike_rates_cs = spike_counts_cs/0.0005;

                spike_counts_us = histcounts(RasterXY_spikes_us(ismember(RasterXY_trials_us, trials)),unit.bins_us);
                spike_counts_us = spike_counts_us/length(trials);
                spike_rates_us = spike_counts_us/0.0005;
                % create filter
                sigma1 = 50; % pick sigma value for the gaussian
                gaussFilter1 = gausswin(6*sigma1 + 1)';
                gaussFilter1 = gaussFilter1 / sum(gaussFilter1); % normalize
                filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter1, 'same');
                filtered_spike_rates_us = conv(spike_rates_us, gaussFilter1, 'same');
                % do an extra gauss filter to smoothen out the lines
                sigma2 = 10; % pick sigma value for the gaussian
                gaussFilter2 = gausswin(6*sigma2 + 1)';
                gaussFilter2 = gaussFilter2 / sum(gaussFilter2); % normalize
                filtered_spike_rates_cs = conv(filtered_spike_rates_cs, gaussFilter2, 'same');
                filtered_spike_rates_us = conv(filtered_spike_rates_us, gaussFilter2, 'same');

                unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates_cs;
                unit.neuron(n).psth_us_reset(c, :) = filtered_spike_rates_us;
            end


            RasterXY_spikes_cs_filtered = unit.neuron(n).RasterXY_cs_sorted_filtered(1,1:3:end);
            RasterXY_trials_cs_filtered = unit.neuron(n).RasterXY_cs_sorted_filtered(2,1:3:end);

            RasterXY_spikes_us_filtered = unit.neuron(n).RasterXY_us_sorted_filtered(1,1:3:end);
            RasterXY_trials_us_filtered = unit.neuron(n).RasterXY_us_sorted_filtered(2,1:3:end);


            for c = 1 : length(unique(stimtimesCS_sorted(:,2)))
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
                spike_counts_cs = histcounts(RasterXY_spikes_cs_filtered(ismember(RasterXY_trials_cs_filtered, trials)),unit.bins_cs);
                spike_counts_cs = spike_counts_cs/sum(ismember(trials, RasterXY_trials_cs_filtered));
                spike_rates_cs = spike_counts_cs/0.0005;

                spike_counts_us = histcounts(RasterXY_spikes_us_filtered(ismember(RasterXY_trials_us_filtered, trials)),unit.bins_us);
                spike_counts_us = spike_counts_us/sum(ismember(trials, RasterXY_trials_cs_filtered));
                spike_rates_us = spike_counts_us/0.0005;
                % create filter
                sigma1 = 50; % pick sigma value for the gaussian
                gaussFilter1 = gausswin(6*sigma1 + 1)';
                gaussFilter1 = gaussFilter1 / sum(gaussFilter1); % normalize
                filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter1, 'same');
                filtered_spike_rates_us = conv(spike_rates_us, gaussFilter1, 'same');
                % do an extra gauss filter to smoothen out the lines
                sigma2 = 10; % pick sigma value for the gaussian
                gaussFilter2 = gausswin(6*sigma2 + 1)';
                gaussFilter2 = gaussFilter2 / sum(gaussFilter2); % normalize
                filtered_spike_rates_cs = conv(filtered_spike_rates_cs, gaussFilter2, 'same');
                filtered_spike_rates_us = conv(filtered_spike_rates_us, gaussFilter2, 'same');

                unit.neuron(n).psth_cs_reset_filtered(c, :) = filtered_spike_rates_cs;
                unit.neuron(n).psth_us_reset_filtered(c, :) = filtered_spike_rates_us;
            end
        end

        path = fileparts(which(FileName_));
        stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'StructEphysData', file_extension = ".mat");
        save(fullfile(path,stamp),'unit')
    end
end

end

function rasterFiltered = rasterFilter(raster, range)
if nargin < 2
    range.min = min(raster(1,:));
    range.max = max(raster(1,:));
end
binWidth = IkUtils.getParams().BinW;
nTrials = max(raster(2,:)) - IkUtils.getParams().raster_spike_height;
binEdges = range.min:binWidth:range.max;
filterMask = zeros(1,length(raster));
for t = 1 : nTrials
    spikeTimesTrial_ = raster(1,raster(2,:) == t);
    rangeMask = ...
        spikeTimesTrial_ >= range.min ...
        & spikeTimesTrial_ <= range.max;
    spikeTimesTrial = spikeTimesTrial_(rangeMask);
    amplitudes = histcounts(spikeTimesTrial, binEdges);
    baseline = mean(amplitudes);
    baserate = baseline/binWidth;
    if baserate > IkUtils.getParams().Sspk_base_filter 
        filterMask(raster(2,:) == t) = 1;
        filterMask(find(raster(2,:) == t) + 1) = 1;
        filterMask(find(raster(2,:) == t) + 2) = 1;
    end
end
rasterFiltered(1,:) = raster(1,find(filterMask));
rasterFiltered(2,:) = raster(2,find(filterMask));
rasterFiltered(3,:) = NaN;
rasterFiltered(3, rasterFiltered(2,:) <= 20.9) = 0;
rasterFiltered(3, rasterFiltered(2,:) > 20.9 & rasterFiltered(2,:) <= 220.9) = 1;
rasterFiltered(3, rasterFiltered(2,:) > 220.9) = 2;
end

function [S, cell_spk, stimtimes] = loadData(mcode, s)

files_stim = recursiveFileSearch('stimtimes*.mat', mcode);
if isempty(files_stim)
    disp("No stimtimes found, this is a Rhymdata mouse or the raw data still needs to be converted.")
    S = [];
    cell_spk = [];
    stimtimes = [];
    return
end

[file_stim, folder_stim] = recursiveFileSearch('stimtimes*.mat', mcode, s);
if isempty(file_stim)
    fprintf("Session %d not found \n", s)
    S = [];
    cell_spk = [];
    stimtimes = [];
    return
elseif numel(folder_stim) > 1
    disp("Too many folders found. Search not specific enough. Check if this is a Rhymdata mouse or the raw data still needs to be converted.")
    S = [];
    cell_spk = [];
    stimtimes = [];
    return
end

stimtimes = load(fullfile(folder_stim, file_stim));
stimtimes = stimtimes.times;

[file_config, folder_config] = recursiveFileSearch("Config_h32_oe*.csv", mcode, s);
if isempty(file_config)
    disp("Config (.csv) file not found")
    S = [];
    cell_spk = [];
    return
end

if numel(folder_config) > 1
    disp("Too many folders found. Search not specific enough. Check if .csv file exists.")
    S = [];
    cell_spk = [];
    return
else
    cell_spk = import_jrc_csv(fullfile(folder_config, file_config));
end

n_neurons = length(cell_spk) - 1;
if n_neurons == 0
    disp("There are no neurons in this session")
    S = [];
    return
end

[file_jrc, folder_jrc] = recursiveFileSearch("Config_h32_oe*jrc.mat", mcode, s);
if numel(folder_jrc) > 1
    disp("Too many folders found. Search not specific enough. Check if jrc.mat file exists.")
    S = [];
    return
else
    S = load(fullfile(folder_jrc, file_jrc));
end

if S.S_clu.nClu ~= n_neurons
    disp("Number of neurons in csv file does not match jrc file")
    S = [];
    return
end

end
%%
trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
spike_counts_cs = histcounts(RasterXY_spikes_cs_filtered(ismember(RasterXY_trials_cs_filtered, trials)),unit.bins_cs);
spike_counts_cs = spike_counts_cs/sum(ismember(trials, RasterXY_trials_cs_filtered));
spike_rates_cs = spike_counts_cs/0.0005;

spike_counts_us = histcounts(RasterXY_spikes_us_filtered(ismember(RasterXY_trials_us_filtered, trials)),unit.bins_us);
spike_counts_us = spike_counts_us/sum(ismember(trials, RasterXY_trials_cs_filtered));
spike_rates_us = spike_counts_us/0.0005;
% create filter
sigma1 = 50; % pick sigma value for the gaussian
gaussFilter1 = gausswin(6*sigma1 + 1)';
gaussFilter1 = gaussFilter1 / sum(gaussFilter1); % normalize
filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter1, 'same');
filtered_spike_rates_us = conv(spike_rates_us, gaussFilter1, 'same');
% do an extra gauss filter to smoothen out the lines
sigma2 = 10; % pick sigma value for the gaussian
gaussFilter2 = gausswin(6*sigma2 + 1)';
gaussFilter2 = gaussFilter2 / sum(gaussFilter2); % normalize
filtered_spike_rates_cs = conv(filtered_spike_rates_cs, gaussFilter2, 'same');
filtered_spike_rates_us = conv(filtered_spike_rates_us, gaussFilter2, 'same');

unit.neuron(n).psth_cs_reset_filtered(c, :) = filtered_spike_rates_cs;
unit.neuron(n).psth_us_reset_filtered(c, :) = filtered_spike_rates_us;
%%
function data = getData(mcode)

P = IkUtils.getParams();
% mouseList = defaultMice();
% dates = mouseList(mouseList.code == mcode).ephysdates;

structFilename = "StructEphysData.mat";
filepattern = "StructEphysData"; fileextension = ".mat";
% structFilePattern = "StructEphysData*.mat";
try
    stamp = nameDateStampFiles(mcode = mcode, s = 1);
    filename_try = fullfile(mcode, P.s(1), filepattern + "_" + stamp + fileextension);
    filename = which(filename_try);
    path_ = fileparts(filename);
    path = erase(path_, P.s(1));
catch
    try
        filename_try = fullfile(mcode, P.s(1), structFilename);
        filename = which(filename_try);
        path_ = fileparts(filename);
        path = erase(path_, P.s(1));
    catch
        disp("Ephys Data Struct not found")
        data = [];
        return
    end
end
mice = defaultMice();
dates = mice([mice.code] == mcode).ephysdates;
if ~isempty(dates)
    numdates = numel(dates);
    sessions = 1:numdates;
else
    disp("No ephys dates found in defaultMice")
    data = [];
    return
end

folders = dir(path); 
folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
folders = folders(ismember(string({folders.name}), P.s));
for s = 1 : length(folders)
    [file_struct, folder_struct] = recursiveFileSearch("StructEphysData*.mat", mcode, s);
    stamp = nameDateStampFiles(mcode = mcode, s = s);
    filename_try = fullfile(mcode, P.s(s), filepattern + "_" + stamp + fileextension);
    filename = which(filename_try);
    if ~exist(filename, 'file')
        filename_try = fullfile(mcode, P.s(s), structFilename);
        filename = which(filename_try);
    end

    if exist(filename, 'file') == 2
        path_ = fileparts(filename);
        path = erase(path_, P.s(1));
    elseif s > 1
        continue
    else
        disp("Ephys Data Struct not found")
        data = [];
        return
    end

    try
        load(filename)
        if ~isempty(unit)            
            data(s) = unit;
        end
    catch
        fprintf("Data file not found for session %d", s)
        data = [];
        if s > 2
            return
        else
            continue
        end
    end
end
if ~exist("data", 'var')
    data = [];
end
end
%%
function reSaveWithStamp(kwargs)
arguments
    kwargs.directory = '/home/i.klinkhamer/Documents/Data/spikeSortingUnits/';
end
directory = kwargs.directory;
P = IkUtils.getParams();

app = getUImain();


while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        % If the UI figure has been closed, exit the loop
        break;
    end
    appData = evalin('base', 'appData');
    mousecodes = appData.mousecodes;
    sessions = appData.sessions;

    % Loop through each file and rename it
    for mcode = mousecodes
        files = dir(fullfile(directory, mcode));
        if ~isempty(files)
%             renameMoveFiles(files, mcode)
        end
        for s = sessions
            files = dir(fullfile(directory, mcode, P.s(s)));
            renameMoveFiles(files, mcode, s)
%             if isempty(files)
%                 continue
%             end
%             stamp = nameDateStampFiles(mcode = mcode, s = s);
% 
%             for i = 1:numel(files)
%                 configFile = 0;
%                 if ~files(i).isdir
%                     % Split filename and extension
%                     [filePath, fileName, fileExt] = fileparts(files(i).name);
% 
%                     parts = split(fileName, '_');
%                     if any(strcmp(parts, "Config"))
%                         configFile = 1;
%                         fileNameNoStamp = erase(files(i).name, stamp);
%                         partsNoStamp_ = split(fileNameNoStamp, '_');
%                         idx = 1;
%                         for p = 1:numel(partsNoStamp_)
%                             if ~isempty(partsNoStamp_{p})
%                                 partsNoStamp{idx} = partsNoStamp_{p};
%                                 idx = idx+1;
%                             end
%                         end
%                         ext = split(fileExt, '.');
%                         partsNoStamp{end} = ext{2};
%                         place = find(~cellfun(@isempty,regexp(partsNoStamp, "dat")) | ~cellfun(@isempty,regexp(partsNoStamp, "bin")));
%                         if isempty(place)
%                             place = find(~cellfun(@isempty,regexp(partsNoStamp, "oe")));
%                         end
%                         part1 = strjoin(partsNoStamp(1:place), '_');
%                         part3_ = string(strjoin(partsNoStamp(place+1:end-1), ""));
%                         part3 = strjoin([part3_, partsNoStamp(end)], ".");
%                         if length(partsNoStamp(place+1:end)) > 1
%                             part3 = "_" + part3;
%                         end
% 
%                         newFileName = part1 + "_" + stamp + part3;
%                         newFileName = fullfile(directory, mcode, P.s(s), newFileName);
%                         if exist(newFileName, 'file') == 2
%                             continue
%                         end
%                     end
%                     if any(strcmp(parts, mcode)) && ~ configFile
%                         continue
%                     end
% 
%                     % Construct the new filename with the new part inserted between original filename and extension
%                     if ~ configFile
%                         newFileName = fullfile(directory, mcode, P.s(s), fileName + "_" + stamp + fileExt);
%                     end
%                     % Move the file to the new filename
%                     movefile(fullfile(directory, mcode, P.s(s), files(i).name), newFileName);
%                 end
%             end
        end
    end
    disp("End of loop. Waiting for new UI input.")
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
            if nargin < 3
                fileNameNoStamp = erase(files(i).name, mcode);
            else
                fileNameNoStamp = erase(files(i).name, nameDateStampFiles(mcode = mcode, s = s));
            end
            partsNoStamp_ = split(fileNameNoStamp, '_');
            idx = 1;
            for p = 1:numel(partsNoStamp_)
                if ~isempty(partsNoStamp_{p})
                    partsNoStamp{idx} = partsNoStamp_{p};
                    idx = idx+1;
                end
            end
            ext = split(fileExt, '.');
            partsNoStamp{end} = ext{2};
            place = find(~cellfun(@isempty,regexp(partsNoStamp, "dat")) | ~cellfun(@isempty,regexp(partsNoStamp, "bin")));
            if isempty(place)
                place = find(~cellfun(@isempty,regexp(partsNoStamp, "oe")));
            end
            file_pattern = strjoin(partsNoStamp(1:place), '_');
            
            file_pattern_end = strjoin(partsNoStamp(place+1:end-1), "");
            file_extension = partsNoStamp(end);            
            if ~isempty(file_pattern_end)%length(partsNoStamp(place+1:end)) > 1
                file_end = strjoin([string(file_pattern_end), file_extension], ".");
                file_end = "_" + file_end;
            else
                file_end = "." + file_extension;               
            end
            newFileNameStamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = file_pattern, file_extension = file_end);
            newFileName = file_pattern + "_" + stamp + file_end;
            newFileName = fullfile(fileFolder, newFileName);
            if exist(newFileName, 'file') == 2
                continue
            end
        end
        if any(strcmp(parts, mcode)) && ~any(strcmp(parts, "Config"))
            continue
        end

        % Construct the new filename with the new part inserted between original filename and extension
        if ~any(strcmp(parts, "Config"))
            newFileName = fullfile(fileFolder, fileName + "_" + stamp + fileExt);
        end
        % Move the file to the new filename
        movefile(fullfile(fileFolder, files(i).name), newFileName);
    end
end

end
%%
function convertData(kwargs)
arguments
    kwargs.previewData = false;
    kwargs.saveConvertedData = true;
    kwargs.saveStimTimestamps = true;
    kwargs.convertAllSessions = true;
    kwargs.convertAllMice = false;
    kwargs.debug = false;
end

mouseCodes = [defaultMice().code];
mouseNames = [defaultMice().name];

ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
mouseNames = mouseNames(ephysMiceMask);
mouseCodes = mouseCodes(ephysMiceMask);

if kwargs.convertAllMice == 0

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to construct a data struct for?")
    fprintf("\n-------------------------------------\n")

    for c = 1 : length(mouseCodes)
        choices(c) = sprintf("%s\t\t%s",mouseNames(c), mouseCodes(c));
    end

    [~, mname_idx] = IkUtils.do_prompt_select_option(choices);
    mcode = mouseCodes(mname_idx);
    mousecodes = mcode;
else
    mousecodes = mouseCodes(kwargs.startMouse:end);
end

if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    mice = defaultMice();
    dates = mice([mice.code] == mcode).ephysdates;
    numdates = sum(~cellfun(@isempty, dates));
    sessions = 1:numdates;
end

for mcode = mousecodes
    for s = sessions
        fprintf("Mouse: %s Session: %d \n", mcode, s);
        path = getPath(mcode, s);

        if kwargs.saveConvertedData == 1
            if ~isempty(dir(fullfile(path, 'RAWdata*.bin')))
                disp("There is already a RAWdata.bin file");
%             elseif ~isempty(dir(fullfile(path, 'continuous*.dat')))
%                 disp("There is already a continuous.dat file");
            else
                saveConvertedData(path, mcode, s);
            end
        end

        if kwargs.previewData == 1
            previewData(path);
        end

        if kwargs.saveStimTimestamps == 1
            if ~isempty(dir(fullfile(path,'stimtimes*.mat'))) && ~isempty(dir(fullfile(path,'timestamps*.mat'))) && kwargs.debug == 0
                disp("There are already timestamp files");
            else
                saveStimTimestamps(path);
            end
        end
    end
end

end

% get path
function path = getPath(mcode, s)
P = IkUtils.getParams();
stamp = nameDateStampFiles(mcode = mcode, s = s);
try
    [~, path] = recursiveFileSearch("all_channels*.events", mcode, s);
catch
    dataFileName = fullfile(mcode, P.s(s), sprintf("all_channels_%s.events", stamp));
    path = fileparts(which(dataFileName));
end
end

% get filepath
function filepath = getFilepath(path, ch)
filepattern = "100_*.continuous";
[files, folders] = recursiveFileSearch(filepattern, mcode, s, parent_directory, max_depth)
try
    fileParts = fileparts(path);
    [~, mcode, num] = fileparts(fileParts);
    mcode = mcode + num;
    reSaveWithStamp(mcode = mcode)
catch
end
filepattern = "100_*.continuous";
dir_contents = dir(fullfile(path, filepattern));
parts_ = split({dir_contents.name}, '_');
if size(parts_,2) > 1
    Parts = reshape(parts_, [size(parts_,2), size(parts_,3)]);
    channelPart = cellfun(@str2double, regexp(Parts(:,2), '\d+', 'match'));
    formatPart_ = regexp(Parts(:,2), '[^\d]+', 'match');
else
    Parts = parts_;
    channelPart = cellfun(@str2double, regexp(Parts(2), '\d+', 'match'));
    formatPart_ = regexp(Parts(2), '[^\d]+', 'match');
end
formatPart = string(cell2mat(formatPart_{1}));
if ~isempty(formatPart) && ch > 16
    ch = ch + 32;
end
filepath = fullfile(dir_contents(channelPart == ch).folder, dir_contents(channelPart == ch).name);
end

% convert data
function saveConvertedData(path, mcode, s)
    filepattern = "100_*.continuous";
    [files, folders] = recursiveFileSearch(filepattern, mcode, s);
    if numel(files) ~= IkUtils.getParams().num_channels
        disp("Number of files found not equal to number of channels on the probe. Returning to main function.")
        return
    end
    if isempty(files)
        disp("Channel.continuous files not found. Returning to main function.")
        return
    end
for i = 1 : IkUtils.getParams().num_channels
    file = files(i);
    folder
    filepath = getFilepath(path, i);
    if isempty(filepath)
        disp("Channel.continuous file not found. Returning to main function.")
        return
    end
    data_ = load_open_ephys_data(filepath);
    if exist('data','var')
        data(i,:) = data_(1:length(data)); %load data
    else
        data(i,:) = data_;
    end
end

try
    [pathmcode, folderS, ~] = fileparts(path);
    [~, mcode_litter, mcode_mouse_end] = fileparts(pathmcode);
    mcode = mcode_litter + mcode_mouse_end;
    stamp = nameDateStampFiles(mcode = mcode, s = str2double(erase(folderS, 'S')));
    fileID = fopen(fullfile(path, sprintf('RAWdata_%s.bin', stamp)),'w');
catch
    fileID = fopen(fullfile(path, 'RAWdata.bin'),'w');
end
fwrite(fileID,data, 'double');
fclose(fileID);
end

% preview data
function previewData(path)
P = getParams();
for ch = 1:P.num_channels
    filepath = getFilepath(path, ch);
    if isempty(filepath)
        continue
    end
    [~, timestamps(ch,:), ~] = load_open_ephys_data(filepath); %load data
end
try
    dir_content = dir(fullfile(path, 'data*.m'));
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

% make timestamps and stimstamps files
function saveStimTimestamps(path)
filepath = getFilepath(path, 1);
if isempty(filepath)
    disp("File not found, returning to main function.")
    return
end
[~,timestamps(1,:), ~] = load_open_ephys_data(filepath); %load data

dir_events = dir(fullfile(path, 'all_channels*.events'));
[~,timestamps2, ~] = load_open_ephys_data(fullfile(path, dir_events.name));

times = timestamps2 - timestamps(1,1);

[pathmcode, folderS, ~] = fileparts(path);
[~, mcode, mcodep2] = fileparts(pathmcode);
mcode = strcat(mcode, mcodep2);
stamp = nameDateStampFiles(mcode = mcode, s = str2double(erase(folderS, 'S')));

save(fullfile(path,sprintf('stimtimes_%s.mat', stamp)),'times')

T = timestamps(:,1);
save(fullfile(path, sprintf('timestamps_%s.mat', stamp)), 'T')
end
%%
function batchProcessTrialsExampleIK(mInput)
close all
debug = 0;
plotfiguresOld = 0;
startMouse = 1;
onlySaveFigs = 1;
% doAll = 1;

P = IkUtils.getParams();
mice = defaultMice();
mouseNames = [mice.name];
mouseCodeList = [mice.code];
if nargin < 1
    mouseCodes = IkUtils.promptMouseNames();
else
    mouseCodes = IkUtils.promptMouseNames(mInput);
end
% if doAll
%     mice = defaultMiceTraining();
%     mouseCodes = [mice.code];
% end

if length(mouseCodes) == 1
    startMouse = 1;
else
    mouseCodes(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT") = [];
end

for m = startMouse:length(mouseCodes)
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
        
        if onlySaveFigs && ~exist(fullfile(figpath, sprintf("EyelidClosure_%s_%s.jpg", mcode, folderName)), 'file') || debug == 1
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
                eyelidpos = upsampleEyelidpos(behavior_trial_data);
                stampeye = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'eyelidpos', file_extension = '.mat');
                                save(fullfile(folder, stampeye), 'eyelidpos');
                if ~exist(fullfile(results_folder, stampeye), 'file')
                    save(fullfile(results_folder, stampeye), 'eyelidpos');
                end
            catch
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
            if ~exist(fullfile(figpath, sprintf("EyelidClosure_%s_%s.jpg", mcode, folderName)), 'file') 
                traces = behavior_trial_data.eyelidpos;
                smoothedTraces = smootheTraces(traces);
                normalized_traces = normalizeEyelidTraces(smoothedTraces, behavior_trial_data);
                fig = plotFigureEyeTraces(behavior_trial_data.tm(1,:), normalized_traces);

                saveas(fig, fullfile(figpath, sprintf("EyelidClosure_%s_%s.jpg", mcode, folderName)));
            end

            if plotfiguresOld
                plotFiguresOld(behavior_trial_data)
            end
        end
    end
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

function eyelidpos = upsampleEyelidpos(behavior_trial_data)
eyelidpos = behavior_trial_data.eyelidpos;
min_val_eye = min(eyelidpos(:));
max_val_eye = max(eyelidpos(:));
eyelidpos = (eyelidpos - min_val_eye) / (max_val_eye - min_val_eye);
% Define the factor by which you want to upsample
upsample_factor = 10;
% Create the original time points
original_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2));
% Create the new time points after upsampling
new_time = linspace(1, size(eyelidpos,2), size(eyelidpos,2)*upsample_factor);
% Interpolate the data

eyelidpos_upsampled = interp1(original_time, eyelidpos.', new_time, 'spline');

smoothedTraces = smootheTraces(eyelidpos_upsampled);
                normalized_traces = normalizeEyelidTraces(smoothedTraces, behavior_trial_data);
min_val_eye = min(eyelidpos_upsampled(:));
max_val_eye = max(eyelidpos_upsampled(:));
eyelidpos = (eyelidpos_upsampled - min_val_eye) / (max_val_eye - min_val_eye);

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
%%

time_points_ = 1:size(eyelidpos, 2);
time_points = [];
for i = 1:size(eyelidpos, 1)
    time_points = vertcat(time_points, time_points_);
end
figure;
plot(time_points, eyelidpos, 'b', 'LineWidth', 0.5);
hold on;
plot(time_points(positive_gradient_mask), eyelidpos(positive_gradient_mask), 'ro', 'MarkerSize', 5);
plot(time_points(negative_gradient_mask), eyelidpos(negative_gradient_mask), 'go', 'MarkerSize', 5);
hold off;
xlabel('Time');
ylabel('Eyelid Trace');
title('Eyelid Trace and Masks');
%     legend('Eyelid Trace', 'Eye Closure', 'Eye Opening', 'Not Closing or Opening Eye', 'Location', 'best');
grid on;
%%

p = IkUtils.getParams();
for m = 1:numel(mouseList)
    mcode = mouseList(m);
    data = getData(mcode);
    for s = 1: length(data)%p.n_sessions
%         [files_stimtimes, folders_stimtimes] = recursiveFileSearch("stimtimes*.mat", mcode, s);
%         stimtimes = load(fullfile(folders_stimtimes, files_stimtimes));
          trial_type_mask = data_eye.type(:,1:p.nTrials) ==1;          
          sessionMask = false(1, size(data_eye.traces,2));
          sessionMask((m-1)*p.nTrials+1:m*p.nTrials) = true;
          mask = trial_type_mask(:)' & sessionMask;
          traces = squeeze(data_eye.traces(s, mask, :));
        for n = 1:length(data(s).neuron)
        
            spikeTimes = data(s).neuron(n).RasterXY_cs(1,1:3:end) * 1000;
            binWidth = IkUtils.getParams().BinW;
            psthRanges = p.psthRanges;
            range = psthRanges.cs_full;
            edges = (range.min:binWidth:range.max) * 1000;
            rangeMask = spikeTimes >= range.min * 1000 & spikeTimes <= range.max * 1000;
            counts = histcounts(spikeTimes(rangeMask), edges);
            figure('Color', 'white','Position', [10 0 1000 1000]);
            hold on
            title(sprintf('%s psth', type))
%             [N1,edges1]=histcounts([data; WTAllData(2).Spikes(1).all_spikes],100);
            a1 = area([-0.2 0.8], [max(counts) max(counts)]);
            a1.FaceColor = [0.5,0.5,0.5];
            a1 = area([0 0.27], [max(counts) max(counts)]);
            a1.FaceColor = [0,0,1];
            a1 = area([0.25 0.27], [max(counts) max(counts)]);
            a1.FaceColor = [0,1,0];
            h = IkUtils.histogramIK('BinCounts', counts, 'BinEdges', edges/1000);
            h.FaceColor = [0.9,0.9,0.9];
            h.FaceAlpha = 0.8;
            plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
            plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
            plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
            plot(data_eye.ts/1000, nanmean(traces)*100, 'LineWidth', 3, 'Color', [0,0,0])
%             plot(linspace(-0.2,0.8,200), data_eye.mean(10,:)*100, 'LineWidth', 3, 'Color', [0,0,0])
            legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
            ylim([min(nanmean(traces)) max(counts)])
            xlim([-0.2 0.8])
            hold off
        end
    end
end
%%

    figure('Color', 'white','Position', [10 0 1000 1000]);
    
    hold on
    title(sprintf('%s psth', type))
[N1,edges1]=histcounts([data; WTAllData(2).Spikes(1).all_spikes],100);
a1 = area([-0.2 0.8], [max(N1) max(N1)]);
a1.FaceColor = [0.5,0.5,0.5];
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1 = area([0.25 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,1,0];
h = histogram('BinCounts', N1, 'BinEdges', edges1);
h.FaceColor = [0.9,0.9,0.9];
h.FaceAlpha = 0.8;
plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
plot(linspace(-0.2,0.8,200), WTmean*1200+100, 'LineWidth', 3, 'Color', [0,0,0])
legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
ylim([0 max(N1)])
xlim([-0.2 2])

MUTpsth = figure('Color', 'white','Position', [10 0 1000 1000]);
hold on
title('MUT psth')
[N1,edges1]=histcounts([MUTAllData(1).Spikes(1).all_spikes; MUTAllData(2).Spikes(1).all_spikes;MUTAllData(3).Spikes(1).all_spikes;MUTAllData(4).Spikes(1).all_spikes;MUTAllData(5).Spikes(1).all_spikes],100);
a1 = area([-0.2 0.8], [max(N1) max(N1)]);
a1.FaceColor = [0.5,0.5,0.5];
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1 = area([0.25 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,1,0];
h = histogram('BinCounts', N1, 'BinEdges', edges1);
h.FaceColor = [0.9,0.9,0.9];
h.FaceAlpha = 0.8;
plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
plot(linspace(-0.2,0.8,200), MUTmean*1000+50, 'LineWidth', 3, 'Color', [0,0,0])
legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
ylim([0 max(N1)])
xlim([-0.2 2])
%%
  trial_type_mask = data_eye.type(:,1:p.nTrials) ==1;          
        sessionMask = false(1, size(data_eye.traces,2));
        sessionMask((m-1)*p.nTrials+1:m*p.nTrials) = true;
        mask = trial_type_mask(:)' & sessionMask;
        traces = squeeze(data_eye.traces(s, mask, :));

        % Create a new figure for each session
        figure('Color', 'white','Position', [10 0 1000 1000]);

        % Plot the mean trace in the upper subplot
        subplot(length(data(s).neuron)+1, 1, 1);
        hold on;
        a1 = area([0 0.27], [max(nanmean(traces))*100 max(nanmean(traces))*100]);
        a1.FaceColor = [0,0,1];a1.FaceAlpha = 0.5;
        a1 = area([0.25 0.27], [max(nanmean(traces))*100 max(nanmean(traces))*100]);
        a1.FaceColor = [0,1,0];
        a1.FaceAlpha = 0.5;
        plot(data_eye.ts/1000, nanmean(traces)*100, 'LineWidth', 3, 'Color', [0,0,0]);
        title(sprintf('Mean Trace - Session %d', s));
        xticklabels({});
        xticks([]);
        xlim([-0.2, 0.75]);
        box('off')
        hold off;

        % Plot the histograms of all neurons in the lower subplot

        for n = 1:length(data(s).neuron)
            subplot(length(data(s).neuron)+1, 1, n+1);
            hold on;
            spikeTimes = data(s).neuron(n).RasterXY_cs(1,1:3:end) * 1000;
            binWidth = IkUtils.getParams().BinW;
            psthRanges = p.psthRanges;
            range = psthRanges.cs_full;
            edges = (range.min:binWidth:range.max) * 1000;
            rangeMask = spikeTimes >= range.min * 1000 & spikeTimes <= range.max * 1000;
            counts = histcounts(spikeTimes(rangeMask), edges);
            a1 = area([0 0.27], [max(counts) max(counts)]);
            a1.FaceColor = [0,0,1];a1.FaceAlpha = 0.5;
            a1 = area([0.25 0.27], [max(counts) max(counts)]);
            a1.FaceColor = [0,1,0];a1.FaceAlpha = 0.5;
            h = IkUtils.histogramIK('BinCounts', counts, 'BinEdges', edges/1000);
            h.FaceColor = [0.9,0.9,0.9];
            h.FaceAlpha = 0.8;
            ylabel('Counts');
            xlim([-0.2, 0.75]);
            box('off')

            if n == length(data(s).neuron) % Only for the bottom subplot
                xlabel('Time (s)');
                xticks(edges(1):0.1:edges(end)/1000); % Set x-axis ticks
            elseif n == 1
                title(sprintf('Histograms of All Neurons - Session %d', s));
            else
                % Remove x-axis ticks and labels for other subplots
                xticklabels({});
                xticks([]);
            end
            hold off;
        end
%%
% IK 28-5-24
% Called by function: saveSummaryFiguresMouseWrapperIK
function plotSummaryFiguresMUTvsWT_IK(mice, gene, loc, outputRoot, loadPrevData)
arguments
    mice
    gene
    loc
    outputRoot
    loadPrevData
end

if ~loadPrevData || ~exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")), 'file')
    % Initialize empty arrays for WT
    allStatsWT = struct;
    allStatsWT.cs = [];
    allStatsWT.us = [];

    % Populate the structures for WT
    for m = 1:length(mice{1})
        mouseCodesGroup = mice{1};
        allStatsWTtemp = histStatsMouseWrapper( ...
            mouseCodesGroup(m), ...
            'event', "cs_facilitation", ...
            'modulating', 1 ...
            );
        if ~isempty(allStatsWTtemp)
            allStatsWT.cs = [allStatsWT.cs allStatsWTtemp.cs];  % Concatenate 'cs' field
            allStatsWT.us = [allStatsWT.us allStatsWTtemp.us];  % Concatenate 'us' field
        end
    end

    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")), "allStatsWT");
else
    load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")));
end

if ~loadPrevData || ~exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")), 'file')
    % Initialize empty arrays for MUT
    allStatsMUT = struct;
    allStatsMUT.cs = [];
    allStatsMUT.us = [];

    % Populate the structures for MUT
    for m = 1:length(mice{2})
        mouseCodesGroup = mice{2};
        allStatsMUTtemp = histStatsMouseWrapper( ...
            mouseCodesGroup(m), ...
            'event', "cs_facilitation", ...
            'modulating', 1 ...
            );
        if ~isempty(allStatsMUTtemp)
            allStatsMUT.cs = [allStatsMUT.cs allStatsMUTtemp.cs];  % Concatenate 'cs' field
            allStatsMUT.us = [allStatsMUT.us allStatsMUTtemp.us];  % Concatenate 'us' field
        end
    end

    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")), "allStatsMUT");
else
    load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")));
end

cs_facilitation_mask_MUT = [allStatsMUT.cs.cs_facilitation] == 1;
cs_facilitation_mask_WT = [allStatsWT.cs.cs_facilitation] == 1;
cs_suppression_mask_MUT = [allStatsMUT.cs.cs_suppression] == 1;
cs_suppression_mask_WT = [allStatsWT.cs.cs_suppression] == 1;
non_mod_mask_cs_MUT = [allStatsMUT.cs.cs_facilitation] == 0 & [allStatsMUT.cs.cs_suppression] == 0;
non_mod_mask_cs_WT = [allStatsWT.cs.cs_facilitation] == 0 & [allStatsWT.cs.cs_suppression] == 0;

us_facilitation_mask_MUT = [allStatsMUT.us.us_facilitation] == 1;
us_facilitation_mask_WT = [allStatsWT.us.us_facilitation] == 1;
us_suppression_mask_MUT = [allStatsMUT.us.us_suppression] == 1;
us_suppression_mask_WT = [allStatsWT.us.us_suppression] == 1;
non_mod_mask_us_MUT = [allStatsMUT.us.us_facilitation] == 0 & [allStatsMUT.us.us_suppression] == 0;
non_mod_mask_us_WT = [allStatsWT.us.us_facilitation] == 0 & [allStatsWT.us.us_suppression] == 0;

% CS

% covariance
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
hold on 
subplot(1, 2, 1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).cov] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).cov] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).cov] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).cov] ...
    [allStatsMUT.cs(cs_suppression_mask_MUT).cov] ...
    [allStatsWT.cs(cs_suppression_mask_WT).cov] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).cov]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).cov]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).cov]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).cov]; ...
    [allStatsMUT.cs(cs_suppression_mask_MUT).cov]; ...
    [allStatsWT.cs(cs_suppression_mask_WT).cov] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('CV')
xlim([min(boxplot_data), max(boxplot_data)])
% xlim([0.002 0.004]);
set(gca,'tickdir','out','box','off');
title(sprintf('CV spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
% covariance
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).cov] ...
    [allStatsWT.us(non_mod_mask_us_WT).cov] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).cov] ...
    [allStatsWT.us(us_facilitation_mask_WT).cov] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).cov]; ...
    [allStatsWT.us(non_mod_mask_us_WT).cov]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).cov]; ...
    [allStatsWT.us(us_facilitation_mask_WT).cov] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('CV')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([0.002 0.004]);
set(gca,'tickdir','out','box','off');
title(sprintf('CV spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallCV_MUTvsWT_%s%s.png', gene, loc)), 'png')
% saveas(gcf, fullfile(outputRoot, sprintf('boxplotallCV_MUTvsWT_%s%s.eps', gene, loc)), 'eps')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallCV_MUTvsWT_%s%s.eps', gene, loc)))


%local covariance
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
subplot(1,2,1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).meanLcov] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).meanLcov] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).meanLcov] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).meanLcov] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).meanLcov]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).meanLcov]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).meanLcov]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).meanLcov] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('LCV')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('LCV spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).meanLcov] ...
    [allStatsWT.us(non_mod_mask_us_WT).meanLcov] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).meanLcov] ...
    [allStatsWT.us(us_facilitation_mask_WT).meanLcov] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).meanLcov]; ...
    [allStatsWT.us(non_mod_mask_us_WT).meanLcov]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).meanLcov]; ...
    [allStatsWT.us(us_facilitation_mask_WT).meanLcov] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('LCV')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('LCV spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallLCV_MUTvsWT_%s%s.png', gene, loc)), 'png')
% saveas(gcf, fullfile(outputRoot, sprintf('boxplotallLCV_MUTvsWT_%s%s.eps', gene, loc)), 'eps')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallLCV_MUTvsWT_%s%s.eps', gene, loc)))

% ISI
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
subplot(1,2,1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).meanISIroi] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).meanISIroi] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).meanISIroi] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).meanISIroi] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).meanISIroi]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).meanISIroi]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).meanISIroi]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).meanISIroi] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('ISI')
xlim([min(boxplot_data), max(boxplot_data)])
%xlim([-0.0001 0.0005]);
set(gca,'tickdir','out','box','off');
title(sprintf('ISI spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).meanISIroi] ...
    [allStatsWT.us(non_mod_mask_us_WT).meanISIroi] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).meanISIroi] ...
    [allStatsWT.us(us_facilitation_mask_WT).meanISIroi] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).meanISIroi]; ...
    [allStatsWT.us(non_mod_mask_us_WT).meanISIroi]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).meanISIroi]; ...
    [allStatsWT.us(us_facilitation_mask_WT).meanISIroi] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('ISI')
xlim([min(boxplot_data), max(boxplot_data)])
% xlim([-0.0001 0.0005]);
set(gca,'tickdir','out','box','off');
title(sprintf('ISI spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallISI_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallISI_MUTvsWT_%s%s.eps', gene, loc)))

% rate
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
subplot(1,2,1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).rate] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).rate] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).rate] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).rate] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).rate]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).rate]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).rate]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).rate] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('freq')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('freq spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).rate] ...
    [allStatsWT.us(non_mod_mask_us_WT).rate] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).rate] ...
    [allStatsWT.us(us_facilitation_mask_WT).rate] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).rate]; ...
    [allStatsWT.us(non_mod_mask_us_WT).rate]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).rate]; ...
    [allStatsWT.us(us_facilitation_mask_WT).rate] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('freq')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('freq spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallfreq_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallfreq_MUTvsWT_%s%s.eps', gene, loc)))


% maxAmpTime
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
subplot(1,2,1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).maxAmpTime] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).maxAmpTime] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).maxAmpTime] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).maxAmpTime] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).maxAmpTime]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).maxAmpTime]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).maxAmpTime]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).maxAmpTime] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('PT')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('PT spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).maxAmpTime] ...
    [allStatsWT.us(non_mod_mask_us_WT).maxAmpTime] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).maxAmpTime] ...
    [allStatsWT.us(us_facilitation_mask_WT).maxAmpTime] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).maxAmpTime]; ...
    [allStatsWT.us(non_mod_mask_us_WT).maxAmpTime]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).maxAmpTime]; ...
    [allStatsWT.us(us_facilitation_mask_WT).maxAmpTime] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('PT')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('PT spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallPT_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallPT_MUTvsWT_%s%s.eps', gene, loc)))

%% SD
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
subplot(1,2,1); hold on
boxplot_data = [ ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).sd] ...
    [allStatsWT.cs(non_mod_mask_cs_WT).sd] ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).sd] ...
    [allStatsWT.cs(cs_facilitation_mask_WT).sd] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.cs(non_mod_mask_cs_MUT).sd]; ...
    [allStatsWT.cs(non_mod_mask_cs_WT).sd]; ...
    [allStatsMUT.cs(cs_facilitation_mask_MUT).sd]; ...
    [allStatsWT.cs(cs_facilitation_mask_WT).sd] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('SD')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('SD spk times CS MUT vs WT'))
hold off

subplot(1,2,2); hold on
boxplot_data = [ ...
    [allStatsMUT.us(non_mod_mask_us_MUT).sd] ...
    [allStatsWT.us(non_mod_mask_us_WT).sd] ...
    [allStatsMUT.us(us_facilitation_mask_MUT).sd] ...
    [allStatsWT.us(us_facilitation_mask_WT).sd] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us(non_mod_mask_us_MUT).sd]; ...
    [allStatsWT.us(non_mod_mask_us_WT).sd]; ...
    [allStatsMUT.us(us_facilitation_mask_MUT).sd]; ...
    [allStatsWT.us(us_facilitation_mask_WT).sd] ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'nonMUT','nonWT','facMUT','facWT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('SD')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('SD spk times US MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotallSD_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotallSD_MUTvsWT_%s%s.eps', gene, loc)))

figurePosition = [100, 100, 1800, 500]; figure('Position', figurePosition);
subplot(1,3,1)
hold on 
boxplot_data = [ ...
    [allStatsMUT.us.rateFull] ...
    [allStatsWT.us.rateFull] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.rateFull]; ...
    [allStatsWT.us.rateFull]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('frequency (Hz)')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('Rate full trial window spk times MUT vs WT'))
hold off

subplot(1,3,2)
hold on 
boxplot_data = [ ...
    [allStatsMUT.us.baserate] ...
    [allStatsWT.us.baserate] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.baserate]; ...
    [allStatsWT.us.baserate]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('frequency (Hz)')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('Baseline spk times MUT vs WT'))
hold off

subplot(1,3,3); hold on
boxplot_data = [ ...
    [allStatsMUT.us.covFull] ...
    [allStatsWT.us.covFull] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.covFull]; ...
    [allStatsWT.us.covFull]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
ylabel('Type')
xlabel('CV')
xlim([min(boxplot_data), max(boxplot_data)])
%     xlim([-0.0001 0.001]);
set(gca,'tickdir','out','box','off');
title(sprintf('Covariance spk times MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotall_rate_and_cov_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotall_rate_and_cov_MUTvsWT_%s%s.eps', gene, loc)))

end
%%
