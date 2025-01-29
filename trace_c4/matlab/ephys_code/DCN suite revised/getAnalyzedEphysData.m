% TODO: Julius: Migrate this code to a central data access library (narainlib)
function data = getAnalyzedEphysData(filelist)

    nSessions = length(filelist);

    for s = 1 : nSessions

        DnDisp(sprintf('Session no : %d',s))
        % Start session loop here....

        % Access the summary file
        % load file >> params, unit and foldername
        DnDisp('Loading file for session')
        load(fullfile(filelist(s).folder,filelist(s).name))

        unit.neuron(1).AveFiring_cs = unit.AveFiring_cs;

        nNeuronsSession(s) = size(unit.neuron(1).AveFiring_cs,1);


        [nConditions, index, ~, ~] = getConditions(unit);
%         [a,b] = GetButter_IK(3); % IK change. I changed a lot of things
%         in this function

        for n = 1:nNeuronsSession(s)

%             idx_cs = find(unit.neuron(n).bins_cs > -0.15 & unit.neuron(n).bins_cs < -0.05);
%             idx_us = find(unit.neuron(n).bins_cs > -0.15 & unit.neuron(n).bins_cs < -0.05); % Same as CS for now

%             base_fr = mean(mean(unit.neuron(n).psth_cs(:,idx_cs)));
            unit.neuron(n).RasterXY_us_sorted = unit.neuron(n).RasterXY_cs;
            unit.neuron(n).RasterXY_us = unit.neuron(n).RasterXY_cs;

            stimtimesON_sorted = unit.stimtimesON_sorted_CS_only;
            RasterXY_spikes_cs = unit.neuron(n).RasterXY_cs(1,1:3:end);
            RasterXY_trials_cs = unit.neuron(n).RasterXY_cs(2,1:3:end);
            RasterXY_spikes_us = unit.neuron(n).RasterXY_us_sorted(1,1:3:end);
            RasterXY_trials_us = unit.neuron(n).RasterXY_us_sorted(2,1:3:end);
           
            for c = 1:nConditions

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
                unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates_cs;    
                unit.neuron(n).psth_us_reset(c, :) = filtered_spike_rates_us;    

%                 unit.neuron(n).psth_intact_cs(c,:) = filter(b,a,unit.neuron(n).psth_cs(c,:));
%                 unit.neuron(n).psth_intact_us(c,:) = filter(b,a,unit.neuron(n).psth_us(c,:));

                % CS- Aligned PSTH --- Mean Reset
%                 fr = filter(b,a,unit.neuron(n).psth_cs(c,:));
%                 unit.neuron(n).psth_cs_reset(c,:) = fr - mean(fr(idx_cs)) + base_fr  ;
% 
%                 fr = filter(b,a,unit.neuron(n).psth_us(c,:));
%                 unit.neuron(n).psth_us_reset(c,:) = fr - mean(fr(idx_us)) + base_fr  ;


                [~,spike_col] = find(unit.neuron(n).binA_cs(c).ba(:,index.complex_spike_range) == 1);
                spike_col = index.complex_spike_range(spike_col);
                spike_times = unit.neuron(n).bins_cs(spike_col);
                unit.neuron(n).CS_st_cond(c).st = spike_times;


            end

            data(s) = unit;

        end

    end
end