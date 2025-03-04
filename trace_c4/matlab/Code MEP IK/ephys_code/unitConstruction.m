%% I.K. 25/10/2023   Taking Spikesorting data outcome and putting it into a nice structure to use.
function unitConstruction(kwargs)
arguments
    kwargs.convertAllSessions = false;
    kwargs.convertAllMice = false;
    kwargs.startMouse = 1;
    kwargs.skipToday = false;
end
P = IkUtils.getParams();
getDateLater = false;
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
    getDateLater = true;
end

for mcode = mousecodes
    if getDateLater
        mice = defaultMice();
        dates = mice([mice.code] == mcode).ephysdates;
        if ~isempty(dates)
            numdates = numel(dates);
            sessions = 1:numdates;
        else
            continue
        end
    end
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mcode, mouseNames(mouseCodes == mcode), s);

        if kwargs.skipToday
            [file_struct, folder_struct] = recursiveFileSearch('StructEphysData*.mat', mcode, s);

            if ~isempty(file_struct)
                dirContentsStruct = dir(fullfile(folder_struct(1), file_struct(1)));
                % Convert modification date string to datetime
                try
                    time_since_modification_date = datetime(dirContentsStruct.date, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss') - datetime('now');
                catch
                    time_since_modification_date = datetime(dirContentsStruct.date, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'Locale', 'nl_NL') - datetime('now');
                end
                currentDate = datetime('today');                
                april30 = datetime(year(currentDate), 4, 30);
                daysDifference = days(april30 - currentDate);
                if days(time_since_modification_date) > daysDifference
                    continue
                end
            end
        end

        [S, cell_spk, stimtimes, path] = loadData(mcode, s);
        if isempty(S) || isempty(cell_spk) || isempty(stimtimes)
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
            stim_times.stimtimes = stimtimes;
            stim_times.lightON = stimtimes(1:2:end);
            stim_times.lightOFF = stimtimes(2:2:end);
            stim_times.CSonlyLightON = stimtimes(21:22:end);
            stim_times.CSonlyLightOFF = stimtimes(22:22:end);

            stimtimesCS = stimtimes(1:2:end); % Take only the odd stimtimes, because those are the times that the light is turning ON.
            stimtimesCS(:,2) = 1;
            stimtimesCS(11:11:end,2) = 0; % IK: CS only trials.

            stim_times.stimtimesCS = stimtimesCS;
            stimtimesCS_sorted = sortrows(stimtimesCS, 2, 'ascend');
            stim_times.stimtimesCS_sorted = stimtimesCS_sorted;

            stimtimes(1:2:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(21:22:end,3) = 1; % IK: CS only trials.
            stimtimes(22:22:end,3) = 1; % IK: CS only trials.

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

            stimtimes_CS_during_US_only_trials = stim_times.USonlyPuffStart - (P.isi / 1000);
            stimtimesCS = sortrows([stim_times.lightON; stimtimes_CS_during_US_only_trials]);


            CS_only_trial_idcs = 12:12:length(stimtimesCS);
            US_only_trial_idcs = 1:12:length(stimtimesCS);
            normal_trial_idcs = setdiff(1:length(stimtimesCS), [CS_only_trial_idcs, US_only_trial_idcs]);
            %             sorted_idcs = [CS_only_trial_idcs, normal_trial_idcs, US_only_trial_idcs];
            %             stimtimesCS_sorted = stimtimesCS(sorted_idcs);

            stimtimesCS(CS_only_trial_idcs, 2) = 0;
            stimtimesCS(normal_trial_idcs, 2) = 1;
            stimtimesCS(US_only_trial_idcs, 2) = 2;

            stim_times.stimtimesCS = stimtimesCS;
            stimtimesCS_sorted = sortrows(stimtimesCS, 2, "ascend");
            stim_times.stimtimesCS_sorted = stimtimesCS_sorted;

            stimtimes(3:4:end,2) = 1; % Times that the light goes on are 1.
            stimtimes(43:44:end,3) = 1; % IK: CS only trials.
            stimtimes(44:44:end,3) = 1; % IK: CS only trials.
            stimtimes(1:44:end, 4) = 1; % IK: US only trials.
            stimtimes(2:44:end, 4) = 1; % IK: US only trials.
        end

        unit.stimtimes = stimtimes;
        unit.stim_times = stim_times;

        n_trials = length(stimtimesCS);

        % extract spiketimes and stimtimes
        for n = 1 : n_neurons
            trial_numbers.spikes(n).trial = zeros(length(unit.neuron(n).spiketimes),1);
            for trial = 1 : n_trials
                trial_numbers.spikes(n).trial(unit.neuron(n).spiketimes > (stimtimesCS_sorted(trial, 1) - P.t_raster_edges) &...
                    unit.neuron(n).spiketimes < (stimtimesCS_sorted(trial, 1) + P.t_raster_edges)) =...
                    trial;
            end
            trial_numbers_neuron = trial_numbers.spikes(n).trial;
            mask_trials = trial_numbers_neuron > 0;
            unit.neuron(n).RasterXY_cs = zeros(2, sum(mask_trials)*3);
            unit.neuron(n).RasterXY_cs(2, 1:3:end) = trial_numbers_neuron(mask_trials);
            unit.neuron(n).RasterXY_cs(2, 2:3:end) = trial_numbers_neuron(mask_trials) + P.raster_spike_height;
            unit.neuron(n).RasterXY_cs(1:2, 3:3:end) = nan;
        end

        for n = 1 : n_neurons
            trial_numbers_neuron = trial_numbers.spikes(n).trial;
            for trial = 1 : n_trials
                unit.neuron(n).trial_spikes_cs(trial).trial = unit.neuron(n).spiketimes(trial_numbers_neuron == trial) - stimtimesCS(trial, 1);
                trial_spikes_cs_sorted = unit.neuron(n).spiketimes(trial_numbers_neuron == trial) - stimtimesCS_sorted(trial, 1);
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial) = trial_spikes_cs_sorted';
                unit.neuron(n).RasterXY_cs(1,unit.neuron(n).RasterXY_cs(2,:) ==  trial + P.raster_spike_height) = trial_spikes_cs_sorted';
            end
            unit.neuron(n).RasterXY_us(1,:) = unit.neuron(n).RasterXY_cs(1,:) + P.t_US_offset;
            unit.neuron(n).RasterXY_us(2,:) = unit.neuron(n).RasterXY_cs(2,:);

            rasters.spikes_cs = unit.neuron(n).RasterXY_cs(1,1:3:end);
            rasters.spikes_us = unit.neuron(n).RasterXY_us(1,1:3:end);
            rasters.trials_cs = unit.neuron(n).RasterXY_cs(2,1:3:end);
            rasters.trials_us = unit.neuron(n).RasterXY_us(2,1:3:end);

            [psth_cs_reset, psth_us_reset] = calcPsthReset(unit, rasters, stimtimesCS_sorted);
            unit.neuron(n).psth_cs_reset = psth_cs_reset;
            unit.neuron(n).psth_us_reset = psth_us_reset;

            unit.neuron(n).RasterXY_cs_filtered = rasterFilter(unit.neuron(n).RasterXY_cs);
            unit.neuron(n).RasterXY_us_filtered = rasterFilter(unit.neuron(n).RasterXY_us);

            rasters.spikes_cs = unit.neuron(n).RasterXY_cs_filtered(1,1:3:end);
            rasters.spikes_us = unit.neuron(n).RasterXY_us_filtered(1,1:3:end);
            rasters.trials_cs = unit.neuron(n).RasterXY_cs_filtered(2,1:3:end);
            rasters.trials_us = unit.neuron(n).RasterXY_us_filtered(2,1:3:end);

            [psth_cs_reset, psth_us_reset] = calcPsthReset(unit, rasters, stimtimesCS_sorted);
            unit.neuron(n).psth_cs_reset_filtered = psth_cs_reset;
            unit.neuron(n).psth_us_reset_filtered = psth_us_reset;

        end

        stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'StructEphysData', file_extension = ".mat");
        save(fullfile(path,stamp),'unit')
    end
end

end


function [S, cell_spk, stimtimes, path] = loadData(mcode, s)

files_stim = recursiveFileSearch('stimtimes*.mat', mcode);
if isempty(files_stim)
    disp("No stimtimes found, this is a Rhymdata mouse or the raw data still needs to be converted.")
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
end

[file_config, folder_config] = recursiveFileSearch("Config_h32_oe*.csv", mcode, s);
if isempty(file_config)
    disp("Config (.csv) file not found")
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
end
if numel(folder_config) > 1
    disp("Too many folders found. Search not specific enough. Check if .csv file exists.")
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
else
    cell_spk = import_jrc_csv(fullfile(folder_config, file_config));
    path = folder_config;
end

n_neurons = length(cell_spk) - 1;

[file_stim, folder_stim] = recursiveFileSearch('stimtimes*.mat', mcode, s);
if isempty(file_stim)
    fprintf("Session %d not found \n", s)
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
elseif numel(folder_stim) > 1
    if ~any(strcmp(folder_stim, folder_config))
        disp("Too many folders found. Search not specific enough. Check if this is a Rhymdata mouse or the raw data still needs to be converted.")
        if ~exist('S', 'var');          S = [];             end
        if ~exist('cell_spk', 'var');   cell_spk = [];      end
        if ~exist('path', 'var');       path = [];          end
        if ~exist('stimtimes', 'var');  stimtimes = [];     end
        return
    else
        file_stim = file_stim(find(folder_stim == folder_config, 1));
        folder_stim = folder_config;
    end
end

stimtimes = load(fullfile(folder_stim, file_stim));
stimtimes = stimtimes.times;

[file_jrc, folder_jrc] = recursiveFileSearch("Config_h32_oe*jrc.mat", mcode, s);

if numel(folder_jrc) > 1
    disp("Too many folders found. Search not specific enough. Check if jrc file exists.")
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
elseif isempty(folder_jrc)
    [file_jrc, folder_jrc] = recursiveFileSearch("Config_h32_oe*.mat", mcode, s);
    try
        file_jrc = file_jrc(1);
        folder_jrc = folder_jrc(1);
        S = load(fullfile(folder_jrc, file_jrc));
    catch
        disp("jrc file not found")
        if ~exist('S', 'var');          S = [];             end
        if ~exist('cell_spk', 'var');   cell_spk = [];      end
        if ~exist('path', 'var');       path = [];          end
        if ~exist('stimtimes', 'var');  stimtimes = [];     end
        return
    end
else
    S = load(fullfile(folder_jrc, file_jrc));
end

if S.S_clu.nClu ~= n_neurons
    disp("Number of neurons in csv file does not match jrc file")
    if ~exist('S', 'var');          S = [];             end
    if ~exist('cell_spk', 'var');   cell_spk = [];      end
    if ~exist('path', 'var');       path = [];          end
    if ~exist('stimtimes', 'var');  stimtimes = [];     end
    return
end

end

function [psth_cs_reset, psth_us_reset] = calcPsthReset(unit, rasters, stimtimes)
P = IkUtils.getParams();
for c = 1 : length(unique(stimtimes(:,2)))
    trials = find(stimtimes(:,2) == P.conditions(c));
    spike_counts_cs = histcounts(rasters.spikes_cs(ismember(rasters.trials_cs, trials)),unit.bins_cs);
    spike_counts_cs = spike_counts_cs/sum(ismember(trials, rasters.trials_cs));
    spike_rates_cs = spike_counts_cs/(unit.bins_cs(1)-unit.bins_cs(2));

    spike_counts_us = histcounts(rasters.spikes_us(ismember(rasters.trials_us, trials)),unit.bins_us);
    spike_counts_us = spike_counts_us/sum(ismember(trials, rasters.trials_cs));
    spike_rates_us = spike_counts_us/(unit.bins_us(1)-unit.bins_us(2));
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

    psth_cs_reset(c, :) = filtered_spike_rates_cs;
    psth_us_reset(c, :) = filtered_spike_rates_us;
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
rasterFiltered(1,:) = raster(1,logical(filterMask));
rasterFiltered(2,:) = raster(2,logical(filterMask));
rasterFiltered(3,:) = NaN;
rasterFiltered(3, rasterFiltered(2,:) <= 20.9) = 0;
rasterFiltered(3, rasterFiltered(2,:) > 20.9 & rasterFiltered(2,:) <= 220.9) = 1;
rasterFiltered(3, rasterFiltered(2,:) > 220.9) = 2;
end

