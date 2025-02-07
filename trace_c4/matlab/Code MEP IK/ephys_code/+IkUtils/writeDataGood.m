function unit = writeDataGood(unit)
P = IkUtils.getParams;
for n = 1:length(unit.neuron)

    %             idx_cs = find(unit.neuron(n).bins_cs > -0.15 & unit.neuron(n).bins_cs < -0.05);
    %             idx_us = find(unit.neuron(n).bins_cs > -0.15 & unit.neuron(n).bins_cs < -0.05); % Same as CS for now

    %             base_fr = mean(mean(unit.neuron(n).psth_cs(:,idx_cs)));

    stimtimesON_sorted = unit.stimtimesON_sorted_CS_only;

    unit.neuron(n).RasterXY_us_sorted = unit.neuron(n).RasterXY_cs; % IK change later
    unit.neuron(n).RasterXY_us = unit.neuron(n).RasterXY_cs;

    RasterXY_spikes_cs = unit.neuron(n).RasterXY_cs(1,1:3:end);
    RasterXY_trials_cs = unit.neuron(n).RasterXY_cs(2,1:3:end);
    RasterXY_spikes_us = unit.neuron(n).RasterXY_us_sorted(1,1:3:end);
    RasterXY_trials_us = unit.neuron(n).RasterXY_us_sorted(2,1:3:end);

    for c = 1:P.n_conditions % IK change it so that it reads conditions from params or something

        trials = find(stimtimesON_sorted(:,2) == P.conditions(c));

        spike_counts_cs = histcounts(RasterXY_spikes_cs(ismember(RasterXY_trials_cs, trials)),unit.bins_cs);
        spike_counts_us = histcounts(RasterXY_spikes_us(ismember(RasterXY_trials_us, trials)),unit.bins_us);

        spike_counts_cs = spike_counts_cs/length(trials);
        spike_counts_us = spike_counts_us/length(trials);
        spike_rates_cs = spike_counts_cs/0.0005;
        spike_rates_us = spike_counts_us/0.0005;

        % create filter
        sigma = 200; % pick sigma value for the gaussian
        gaussFilter = gausswin(6*sigma + 1)';
        gaussFilter = gaussFilter / sum(gaussFilter); % normalize
        filtered_spike_rates_cs = conv(spike_rates_cs, gaussFilter, 'same');
        filtered_spike_rates_us = conv(spike_rates_us, gaussFilter, 'same');
        unit.neuron(n).psth_cs_reset = [];
        unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates_cs;
        unit.neuron(n).psth_us_reset(c, :) = filtered_spike_rates_us;



    end

end



end