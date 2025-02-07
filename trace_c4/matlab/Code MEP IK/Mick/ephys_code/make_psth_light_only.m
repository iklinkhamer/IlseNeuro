

%path ='/home/mick/Desktop/0112507_2023-05-11_08-38-54_1/RecordNode105/experiment1/recording1/continuous/Rhythm_FPGA-100.0/';
path = '/home/mick/Desktop/024404/01/';
t_pre = 500;
t_post = 250; %for CR detection
dur = 2800;
t_psth = 1000;
bin_psth = 150; 
bin_corr = 250;

%cell_spk = import_jrc_csv(strcat(path,'RAWdata_DBC_3.1-64-H2.csv'));
cell_spk = import_jrc_csv(strcat(path,'Config_h32_oe.csv'));

%% psth
load(strcat(path,'stimtimes.mat'));
%times = times(times<=3300);
times = [times zeros(length(times),1)];

times =  times(1:11:end,:);

times(1:2:end,2) = 1;

% extract spiketimes and stimtimes
for j = 2:size(cell_spk,2)
   cell_spk(j).spikes = zeros(length(cell_spk(j).t),1);
for i=1:2:length(times)
    cell_spk(j).spikes(cell_spk(j).t>(times(i)-0.500) & cell_spk(j).t<(times(i)+2.00))=(i+1)/2;
end
end

spikes=nan(length(times)/2,100000,size(cell_spk,2)-1);
for j=1:(size(cell_spk,2)-1)
    jj=2:size(cell_spk,2);
for i=1:length(times)/2
    spikes(i,1:length(cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)),j)=cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)-times((2*i-1));
end
end

psth = figure('Color', 'white','Position', [10 0 1500 700]);
for i=1:(size(cell_spk,2)-1)
subplot(5,5,i)
hold on
title(i)
[N1,edges1] = histcounts(spikes(:,:,i),100);
h = histogram('BinCounts', N1, 'BinEdges', edges1, 'FaceColor', [0.9, 0.9, 0.9], 'FaceAlpha', 0.8);

% Add presentation markers for CS and US
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1.FaceAlpha = 0.3;



% Plot lines to mark the presentation of CS and US
hold on
plot([0 0.27], [max(N1)*1.1 max(N1)*1.1], 'LineWidth', 2, 'Color', [0,0,1])




% Set axis limits and labels
hold off
xlim([-0.5 2])
ylim([0 max(N1)*1.2])
xlabel('Time (s)', 'FontSize', 14)
ylabel('Spike Count', 'FontSize', 14)




% Add a title to the plot
title(sprintf(' CS of neuron %d', i), 'FontSize', 16)
if i ==1
legend( 'CS presented', 'US presented', 'Location', 'southeast')
end
end



% Save psth and spikes_count in file for CS_present 
Spikes(1).all_spikes = nan((size(cell_spk,2)-1)*20,100000);

for i=1:j
    Spikes(i).per_cell=spikes(:,:,i);
    Spikes(1).all_spikes(20*i-19:20*i,1:100000)=spikes(:,:,i);
end

save(strcat(path,'spikes_CS.mat'),'Spikes')
print(psth, fullfile(path, 'psth-percell_CS.png'), '-dpng')
