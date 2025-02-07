%% I.K. 8/9/23   PSTH for combined US/CS
function makeRasterPlots(kwargs)
arguments
    kwargs.loadDataFromDATA = true;
    kwargs.plotAllSessions = true;
    kwargs.plotAllMice = false;
    kwargs.sorted_CS_only_trials = true;
    kwargs.startMouse = 1;
end
close all
addpath('/media/mick/DATA/Ilse/convertedData/')

P = getParams();

if kwargs.plotAllMice == 0
    prompt = "Enter mouse index: ";
    mousenames = string(input(prompt, 's'));
    if any(contains(mousenames, P.mouseNames))
        mousenames = P.mouseList(find(P.mouseNames == mousenames));
    end
else
    mousenames = P.mouseList;
    mousenames = mousenames(kwargs.startMouse:end);
end

if kwargs.plotAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    sessions = 1:P.n_sessions;
end

for mname = mousenames
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mname, P.mouseNames(find(P.mouseList == mname)), s);
        unit = getData(mname, s);
        if isempty(unit)
            continue
        end
        
        n_neurons = length(unit.neuron);

        % psth figure

        for n = 1 : n_neurons
            fig = figure(1);
            set(fig, 'Color', 'white','Position', [0 2000 1800 1000]);
            subplot(1,3,1)
            hold on

            [Spike_counts,bins] = histcounts(unit.neuron(n).RasterXY_cs(1, 1:3:end), 100);


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
            title(sprintf('Mouse: %s, %s  PSTH of neuron %d', mname, P.mouseNames(find(P.mouseList == mname)), n), 'FontSize', 16)
            if n ==1
                legend( 'CS presented', 'US presented', 'Location', 'southeast')
            end


            % CS- Aligned Raster
            %         raster = figure('Color', 'white','Position', [10 0 1500 700]);
            subplot(1,3,2)
            hold on

            title(sprintf('CS Aligned: Sess %d  N  %d', s, unit.irc_trueID(n)))
%             RasterXY_trials = unit.neuron(n).RasterXY_cs(2,:);
%             RasterXY_trials_CS_only_first = unit.sorted_CS_only_trials(unit.neuron(n).RasterXY_cs(2, 1:3:end));
%             RasterXY_trials(1:3:end) = RasterXY_trials_CS_only_first;
%             RasterXY_trials(2:3:end) = RasterXY_trials_CS_only_first + 0.9;
%             unit.neuron(n).RasterXY_cs = [unit.neuron(n).RasterXY_cs(1,:); RasterXY_trials];
            b1 = area([-0.5 2], [20.9 20.9]);
            b1.FaceColor = [0,0.5,1];
            b1.FaceAlpha = 0.3;
            b1.EdgeColor = [0, 0.5, 1];
            b1.EdgeAlpha = 0.3;
            if kwargs.sorted_CS_only_trials == 0
                plot(unit.neuron(n).RasterXY_cs(1,:),unit.neuron(n).RasterXY_cs(2,:),'Color',[0.2 0.2 0.2])
            elseif kwargs.sorted_CS_only_trials == 1
                plot(unit.neuron(n).RasterXY_cs(1,:),unit.neuron(n).RasterXY_cs(2, :),'Color',[0.2 0.2 0.2])
            end
            mn = 1; mx = P.n_trials + 0.9;
            plot([0 0], [mn mx],'k--')

            plot([P.t_US P.t_US], [mn mx],'k--')

            xlim([-0.5 2]); ylim([1 mx]);
          
            hold off

            % CS-aligned lines
            subplot(1,3,3), hold on
            
            stimtimesCS_sorted = unit.stimtimesCS;
            RasterXY_spikes = unit.neuron(n).RasterXY_cs(1,1:3:end);
            RasterXY_trials = unit.neuron(n).RasterXY_cs(2,1:3:end);
            idx_cs = find(unit.neuron(n).RasterXY_cs(1,1:3:end) > -0.15 & unit.neuron(n).RasterXY_cs(1,1:3:end) < -0.05);
            base_fr = mean(mean(RasterXY_spikes(1,idx_cs)));
            for c = 1 : P.n_conditions
                trials = find(stimtimesCS_sorted(:,2) == P.conditions(c));
                for trial = 1 : length(trials)
                    [Spike_counts,bins] = histcounts(RasterXY_spikes(ismember(RasterXY_trials, trials)),100);
                end
                mask_cs_bins = bins > -0.15 & bins < -0.05;
                [a, b] = GetButter();
                fr = filter(b,a,Spike_counts);
                unit.neuron(n).lines_cs_reset(c).fr = fr - mean(fr(mask_cs_bins)) + base_fr  ;

                Spike_counts = Spike_counts./length(trials);
                % create filter
                sigma = 2; % pick sigma value for the gaussian
                gaussFilter = gausswin(6*sigma + 1)';
                gaussFilter = gaussFilter / sum(gaussFilter); % normalize
                filtered_spike_counts = conv(Spike_counts, gaussFilter, 'same');


                % plot(bins(1:end-1), Spike_counts)
                plot(bins(1:end-1), filtered_spike_counts)
                ymax(c) = max(Spike_counts)*1.2;
            end
            ymax = max(ymax);
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

% save(strcat(config_path,'spikesCSUS.mat'),'Spikes')
%print(psth, fullfile(config_path, 'psth.png'), '-dpng')


end

function unit = getData(mname, s)
P = getParams();
directory = P.dirHome;
config_path = fullfile(directory, mname, P.s(s));

if ~isfolder(config_path) && s == 2
    unit = [];
    disp("Only one session found")
    return
end

if ~isempty(dir(fullfile(config_path, "StructEphysData.mat")))
    epdir = dir(fullfile(config_path, "StructEphysData.mat"));
    filename = fullfile(epdir.folder, epdir.name);
    load(filename)
else
    unit = [];   
    disp("Data file not found")
    return
end

end
