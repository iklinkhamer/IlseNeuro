function data = getData(mcode, n_sessions, s_start, calc_psth_reset)
arguments
    mcode 
    n_sessions = getParams().n_sessions;    
    s_start = 1;
    calc_psth_reset = false;
end

P = getParams();
directory = P.pathSpikeSortingHome;
s_end = s_start + n_sessions - 1;
sessions  = s_start : s_end;
for s = 1 : n_sessions
    s2 = sessions(s);
    config_path = fullfile(directory, mcode, P.s(s2));

    if ~isfolder(config_path) && s2 == 2
        disp("Only one session found")
        data = [];
        return
    end

    if ~isempty(dir(fullfile(config_path, "StructEphysData.mat")))
        epdir = dir(fullfile(config_path, "StructEphysData.mat"));
        filename = fullfile(epdir.folder, epdir.name);
        load(filename)

        if calc_psth_reset
            for n = 1:length(unit.neuron)

                stimtimesON_sorted = unit.stimtimesON_sorted_CS_only;
                RasterXY_spikes_cs = unit.neuron(n).RasterXY_cs(1,1:3:end);
                RasterXY_trials_cs = unit.neuron(n).RasterXY_cs(2,1:3:end);
                RasterXY_spikes_us = unit.neuron(n).RasterXY_us_sorted(1,1:3:end);
                RasterXY_trials_us = unit.neuron(n).RasterXY_us_sorted(2,1:3:end);

                for c = 1:P.n_conditions

                    trials = find(stimtimesON_sorted(:,2) == P.conditions(c));

                    spike_counts_cs = histcounts(RasterXY_spikes_cs(ismember(RasterXY_trials_cs, trials)),unit.bins_cs);
                    spike_counts_us = histcounts(RasterXY_spikes_us(ismember(RasterXY_trials_us, trials)),unit.bins_us);

                    spike_counts_cs = spike_counts_cs/length(trials);
                    spike_counts_us = spike_counts_us/length(trials);
                    spike_rates_cs = spike_counts_cs/0.0005;
                    spike_rates_us = spike_counts_us/0.0005;

                    % create filter
                    sigma = 100; % pick sigma value for the gaussian
                    gaussFilter = gausswin(6*sigma + 1)';
                    gaussFilter = gaussFilter / sum(gaussFilter); % normalize
                    filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter, 'same');
                    filtered_spike_rates_us = conv(spike_rates_us, gaussFilter, 'same');
                    unit.neuron(n).psth_cs_reset(c, :) = [];
                    unit.neuron(n).psth_us_reset(c, :) = [];
                    unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates_cs;
                    unit.neuron(n).psth_us_reset(c, :) = filtered_spike_rates_us;
                end
            end
        end
        data(s) = unit;
    else
        fprintf("Data file not found for session %d", s2)
        data = [];
        return
    end
end
if ~exist("data", "var")
    data = [];
end

end