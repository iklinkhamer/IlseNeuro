clear all
%PSTH for combined US/CS
%%
p = getParams();
%path = '/home/mick/Desktop/024404/01/';
csv_path = fullfile(p.pathSpikeSortingHome, "20-MI19442-03", "S1");

        if ~isempty(dir(fullfile(csv_path, "*.csv")))
            epdir = dir(fullfile(csv_path, "*.csv"));
            csv_filename = fullfile(epdir.folder, epdir.name);
        else
            disp("csv file not found")
            return
        end
cell_spk = import_jrc_csv(csv_filename);

data_path = fullfile(p.pathEphysDataDATA, "20-MI19442-03", "S1");
%% psth and save US and CS
load(fullfile(data_path,'stimtimes.mat'));
times = times(times<=3500);
%Times = 220
Times = [times zeros(length(times),1)]; 
times(1:11:end,:) = [];
times(1:2:end,2) = 1;


for n = 2:size(cell_spk,2)
   cell_spk(n).spikes = zeros(length(cell_spk(n).t),1);
for trial=1:2:length(times)
    cell_spk(n).spikes(cell_spk(n).t>(times(trial)-0.500) & cell_spk(n).t<(times(trial)+2.00))=(trial+1)/2;
end
end

spikes=nan(length(times)/2,10000,size(cell_spk,2)-1);
n_neurons = (size(cell_spk,2)-1);
n_trials = length(times)/2;
for n=1:n_neurons
    n2=2:size(cell_spk,2);
for trial=1:n_trials
    spikes(trial,1:length(cell_spk(n2(n)).t(cell_spk(n2(n)).spikes==trial)),n)=cell_spk(n2(n)).t(cell_spk(n2(n)).spikes==trial)-times((2*trial-1));
end
end

psth = figure('Color', 'white','Position', [10 0 1500 700]);
for n=1:(size(cell_spk,2)-1)
subplot(5,5,n)
hold on
title(n)
[N1,edges1] = histcounts(spikes(:,:,n)/200,100);
h = IkUtils.histogramIK('BinCounts', N1, 'BinEdges', edges1, 'FaceColor', [0.9, 0.9, 0.9], 'FaceAlpha', 0.8);

% Add presentation markers for CS and US
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1.FaceAlpha = 0.3;

a2 = area([0.25 0.27], [max(N1) max(N1)]);
a2.FaceColor = [0,1,0];
a2.FaceAlpha = 0.3;

% Plot lines to mark the presentation of CS and US
hold on
plot([0 0.27], [max(N1)*1.1 max(N1)*1.1], 'LineWidth', 2, 'Color', [0,0,1])
plot([0.25 0.27], [max(N1)*1.1 max(N1)*1.1], 'LineWidth', 2, 'Color', [0,1,0])



% Set axis limits and labels
hold off
xlim([-0.5 2])
ylim([0 max(N1)*1.2])
xlabel('Time (s)', 'FontSize', 14)
ylabel('Spike Count', 'FontSize', 14)



% Add a title to the plot
title(sprintf('PSTH of neuron %d', n), 'FontSize', 16)
if n ==1
legend( 'CS presented', 'US presented', 'Location', 'southeast')
end
end


%% saving
Spikes(1).all_spikes = nan((size(cell_spk,2)-1)*200,10000);
for trial=1:n_neurons
    Spikes(trial).per_cell=spikes(:,:,trial);
    Spikes(1).all_spikes(200*trial-199:200*trial,1:10000)=spikes(:,:,trial);
end

save(strcat(csv_path,'spikescsUS.mat'),'Spikes')
print(psth, fullfile(csv_path, 'psth.png'), '-dpng')

