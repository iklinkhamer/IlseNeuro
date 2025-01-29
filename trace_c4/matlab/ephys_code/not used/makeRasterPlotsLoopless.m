%% I.K. 8/9/23   PSTH for combined US/CS
function makeRasterPlotsLoopless(mousenames, sessions)
arguments
    mousenames = mouseList();
    sessions = 1:2;
end
%close all
P = IkUtils.getParams();

if any(contains(mousenames, P.mouseNames))
    for i = 1:length(mousenames)
        name = mousenames(i);
        mousenames(i) = P.mouseList(find(P.mouseNames == name));
    end
end

for mname = mousenames
    data = getData(mname);
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mname, P.mouseNames(find(P.mouseList == mname)), s);
        %unit = getData(mname, s);
        try
            unit = data(s);
            unit.bins_cs = -(P.t_pre_trial-0.0005):0.0005:(P.t_post_trial-0.0005);
        catch
            continue
        end
        if isempty(unit)
            continue
        end

        n_neurons = length(unit.neuron);

        % psth figure

        for n = 1 : n_neurons
            %ax = axes('Position', [0 2000 1800 1000]);
            %fig = figure(1);
            fig = figure('Position', [0 2000 1800 1000]);
            set(fig, 'Color', 'white');
            subplot(1,3,1, Parent = fig);
            hold on

            [spike_counts,bins] = histcounts(unit.neuron(n).RasterXY_cs(1, 1:3:end), 100);


            h = IkUtils.histogramIK('BinCounts', spike_counts, 'BinEdges', bins, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0.9, 0.9, 0.9], 'FaceAlpha', 0.8);

            % Add presentation markers for CS and US
            a1 = area([0 0.27], [max(spike_counts) max(spike_counts)]);
            a1.FaceColor = [0,0,1];
            a1.FaceAlpha = 0.3;

            a2 = area([0.25 0.27], [max(spike_counts) max(spike_counts)]);
            a2.FaceColor = [0,1,0];
            a2.FaceAlpha = 0.3;

            % Plot lines to mark the presentation of CS and US
            hold on
            plot([0 0.27], [max(spike_counts)*1.1 max(spike_counts)*1.1], 'LineWidth', 2, 'Color', [0,0,1])
            plot([0.25 0.27], [max(spike_counts)*1.1 max(spike_counts)*1.1], 'LineWidth', 2, 'Color', [0,1,0])

            ymax = max(spike_counts)*1.2;

            % Set axis limits and labels
            hold off
            xlim([-0.5 2])
            ylim([0 ymax])
            xlabel('Time (s)', 'FontSize', 14)
            ylabel('Spike Count', 'FontSize', 14)



            % Add a title to the plot
            title(sprintf('Mouse: %s, %s  PSTH of neuron %d', mname, P.mouseNames(find(P.mouseList == mname)), n), 'FontSize', 16)
            if n ==1
                legend( 'CS presented', 'US presented', 'Location', 'southeast')
            end


            % CS- Aligned Raster
            subplot(1,3,2)
            hold on

            title(sprintf('CS Aligned: Sess %d  N  %d', s, unit.irc_trueID(n)))
            b1 = area([-0.5 2], [20.9 20.9]);
            b1.FaceColor = [0,0.5,1];
            b1.FaceAlpha = 0.3;
            b1.EdgeColor = [0, 0.5, 1];
            b1.EdgeAlpha = 0.3;

            plot(unit.neuron(n).RasterXY_cs(1,:),unit.neuron(n).RasterXY_cs(2, :),'Color',[0.2 0.2 0.2])

            mn = 1; mx = P.n_trials + 0.9;
            plot([0 0], [mn mx],'k--')

            plot([P.t_US P.t_US], [mn mx],'k--')

            xlim([-0.5 2]); ylim([1 mx]);

            hold off

            % CS-aligned lines
            subplot(1,3,3), hold on
            bins = unit.bins_cs;
            try
                stimtimesCS_sorted = unit.stimtimesCS;
            catch
                stimtimesCS_sorted = unit.stimtimesON_sorted_CS_only;
            end
            RasterXY_spikes = unit.neuron(n).RasterXY_cs(1,1:3:end);
            RasterXY_trials = unit.neuron(n).RasterXY_cs(2,1:3:end);

            for c = 1 : P.n_conditions
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));

                spike_counts = histcounts(RasterXY_spikes(ismember(RasterXY_trials, trials)),unit.bins_cs);

                spike_counts = spike_counts./length(trials);
                spike_rates = spike_counts/0.0005;

                % create filter

                sigma = 125; % pick sigma value for the gaussian

                gaussFilter = gausswin(6*sigma + 1)';
                gaussFilter = gaussFilter / sum(gaussFilter); % normalize
                filtered_spike_rates = conv(spike_rates, gaussFilter, 'same');
               
                % filter again to smooth out the lines
                sigma = 20; % pick sigma value for the gaussian
                gaussFilter = gausswin(6*sigma + 1)';
                gaussFilter = gaussFilter / sum(gaussFilter); % normalize
                filtered_spike_rates = conv(filtered_spike_rates, gaussFilter, 'same');
                

                plot(bins(1:end-1), filtered_spike_rates)

                unit.neuron(n).psth_cs_reset = [];
                unit.neuron(n).psth_cs_reset(c, :) = filtered_spike_rates;

            end
            ymax = max([100, max(filtered_spike_rates)+10]);
            mn = 0; mx = ymax;
            plot([0 0], [mn mx],'k--')

            plot([P.t_US P.t_US], [mn mx],'k--')

            xlim([-0.5 2]); ylim([0 ymax]);
            legend("CS only", "CS-US paired");
            hold off


            pause

            figure(1),subplot(1,2,1), cla; figure(1),subplot(1,2,2), cla;
        end
    end
end


%% saving

end
% 
% function unit = getData(mname, s)
% P = getParams2();
% directory = P.dirHome;
% config_path = fullfile(directory, mname, P.s(s));
% 
% if ~isfolder(config_path) && s == 2
%     unit = [];
%     disp("Only one session found")
%     return
% end
% 
% if ~isempty(dir(fullfile(config_path, "StructEphysData.mat")))
%     epdir = dir(fullfile(config_path, "StructEphysData.mat"));
%     filename = fullfile(epdir.folder, epdir.name);
%     load(filename)
% else
%     unit = [];   
%     disp("Data file not found")
%     return
% end
% 
% end
