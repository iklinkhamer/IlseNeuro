%% Ilse 8/9/23   PSTH for combined US/CS
function psthUSCS(kwargs)
arguments
    kwargs.n_sessions = 2;
    kwargs.loadDataFromDATA = true;
end
addpath('/media/mick/DATA/Ilse/convertedData/')

prompt = "Enter mouse index: ";
mname = input(prompt, 's');
prompt2 = "Enter session number: ";
session = input(prompt2);

if session == 1
    session_n = "S1";
elseif session == 2
    session_n = "S2";
end

P = getParams();
directory = P.dirHome;
path_end = mname + "/" + session_n + "/";
csv_path = directory + path_end;

if kwargs.loadDataFromDATA == 1
    dataPath = P.dataPathDATA + path_end;
else
    dataPath = P.dataPathHome + path_end;
end

filename = getFilename(csv_path);
cell_spk = import_jrc_csv(strcat(filename));

%% psth and save US and CS
stimtimes = load(strcat(dataPath,'stimtimes.mat'));
stimtimes = stimtimes.times(stimtimes.times <= 3500);
%Times = 220
% Times = [stimtimes zeros(length(stimtimes),1)]; 
stimtimes(1:11:end,:) = [];
stimtimes(1:2:end,2) = 1;


for j = 2:size(cell_spk,2)
    cell_spk(j).spikes = zeros(length(cell_spk(j).t),1);
    for i=1:2:length(stimtimes)
        cell_spk(j).spikes(cell_spk(j).t>(stimtimes(i)-0.500) & cell_spk(j).t<(stimtimes(i)+2.00))=(i+1)/2;
    end
end

spikes=nan(length(stimtimes)/2,10000,size(cell_spk,2)-1);
for j=1:(size(cell_spk,2)-1)
    jj=2:size(cell_spk,2);
    for i=1:length(stimtimes)/2
        spikes(i,1:length(cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)),j)=cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)-stimtimes((2*i-1));
    end
end

psth = figure('Color', 'white','Position', [10 0 1500 700]);
for i=1:(size(cell_spk,2)-1)
    subplot(5,5,i)
    hold on
    title(i)
    [N1,edges1] = histcounts(spikes(:,:,i)/200,100);
    h = histogram('BinCounts', N1, 'BinEdges', edges1, 'FaceColor', [0.9, 0.9, 0.9], 'FaceAlpha', 0.8);

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

    ymax = max(N1)*1.2;

    % Set axis limits and labels
    hold off
    xlim([-0.5 2])
    ylim([0 ymax])
    xlabel('Time (s)', 'FontSize', 14)
    ylabel('Spike Count', 'FontSize', 14)



    % Add a title to the plot
    title(sprintf('PSTH of neuron %d', i), 'FontSize', 16)
    if i ==1
        legend( 'CS presented', 'US presented', 'Location', 'southeast')
    end
end


%% saving
Spikes(1).all_spikes = nan((size(cell_spk,2)-1)*200,10000);
for i=1:j
    Spikes(i).per_cell=spikes(:,:,i);
    Spikes(1).all_spikes(200*i-199:200*i,1:10000)=spikes(:,:,i);
end

save(strcat(csv_path,'spikescsUS.mat'),'Spikes')
print(psth, fullfile(csv_path, 'psth.png'), '-dpng')


end

function filename = getFilename(csv_path)
if isfile(csv_path + "Config_h32_oe_bin.csv")
    filename = csv_path + "Config_h32_oe_bin.csv";
elseif isfile(csv_path + "Config_h32_oe_dat.csv")
    filename = csv_path + "Config_h32_oe_dat.csv";
elseif isfile(csv_path + "Config_h32_oe.csv")
    filename = csv_path + "Config_h32_oe.csv";
else
    disp("csv file not found")
    return
end
end
