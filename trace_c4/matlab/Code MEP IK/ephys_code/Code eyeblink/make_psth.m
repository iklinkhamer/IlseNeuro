clear all
close all

%%
path = 'D:\Nynke\Ephys-data\MI19442-06F_2021-05-31_13-28-33_2\Record Node 101\';

t_pre = 500;
t_post = 250; %for CR detection
dur = 2000;
t_psth = 1000;
bin_psth = 150; 
bin_corr = 250;

cell_spk = import_jrc_csv(strcat(path,'RAWdata_cnt_edge_32.csv'));

load('D:\Nynke\neuroblink-data\MI1944206\trialdata.mat')
day = 10;
%% Firing freq, ISI, CV, CV2
t = 1800; %seconds
for i = 2:10 %always skip 1, because cluster 0 is noise channel
FF(i) = numel(cell_spk(i).t)/t; 
ISI = cell_spk(i).t(2:end)-cell_spk(i).t(1:(end-1));
CV(i) = std(ISI)/mean(ISI); % std(ISI)/mean(ISI)
CV2(i)= mean(2*abs(ISI(2:end)-ISI(1:(end-1)))./(ISI(2:end)+ISI(1:(end-1)))); % 2*abs(ISI(2)-ISI(1))/(ISI(2)+ISI(1))
end
%% psth
load(strcat(path,'stimtimes.mat'));
times = times(times<=1800);
times = [times zeros(length(times),1)];
times(1:2:end,2) = 1;


for j = 2:size(cell_spk,2)
   cell_spk(j).spikes = zeros(length(cell_spk(j).t),1);
for i=1:2:length(times)
    cell_spk(j).spikes(cell_spk(j).t>(times(i)-0.200) & cell_spk(j).t<(times(i)+2.00))=(i+1)/2;
end
end

spikes=nan(length(times)/2,10000,size(cell_spk,2)-1);
for j=1:(size(cell_spk,2)-1)
    jj=2:size(cell_spk,2);
for i=1:length(times)/2
    spikes(i,1:length(cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)),j)=cell_spk(jj(j)).t(cell_spk(jj(j)).spikes==i)-times((2*i-1));
end
end

%%
psth = figure('Color', 'white','Position', [10 0 1500 700]);
for i=1:(size(cell_spk,2)-1)
subplot(2,5,i)
hold on
title(i)
[N1,edges1]=histcounts(spikes(:,:,i),40);
a1 = area([-0.2 0.8], [max(N1) max(N1)]);
a1.FaceColor = [0.5,0.5,0.5];
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1 = area([0.25 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,1,0];
[N1,edges1]=histcounts(spikes(:,:,i),50);
h = histogram('BinCounts', N1, 'BinEdges', edges1);
h.FaceColor = [0.9,0.9,0.9];
h.FaceAlpha = 0.8;
plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
if i ==1
%legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
end
end

% subplot(1,2,2)
% hold on
% a1 = area([-0.2 0.8], [350 350],'LineStyle',':');
% a1.FaceColor = [0.5,0.5,0.5];
% a1 = area([0 0.27], [350 350],'LineStyle',':');
% a1.FaceColor = [0,0,1];
% a1 = area([0.25 0.27], [350 350],'LineStyle',':');
% a1.FaceColor = [0,1,0];
% [N2,edges2]=histcounts(spikes(:,:,2),40);
% h2 = histogram('BinCounts', N2, 'BinEdges', edges2);
% h2.FaceColor = [0.9,0.9,0.9];
% h2.FaceAlpha = 0.8;

%% correlation
%load('E:\Nynke\neuroblink-data\AllTrials.mat')

% baseline = mean(N1(1:8));
% N1_corrected = N1/baseline;
% 
% figure
% subplot(2,2,1)
% %plot(-AllTrials.MUTmean(10,round(linspace(1,200,40))),linspace(-200,800,40));
% ylim([-200 800])
% subplot(2,2,4)
% plot(linspace(-200,800,40), N1_corrected)
% 
% 
% %R=corrcoef([N1_corrected AllTrials.MUTmean(10,round(linspace(1,200,40)))]);
% 
% for j=1:2
%     jj=2:size(cell_spk,2);
% for i=1:2:length(times)
%    spikerateCS(j,(i+1)/2)=sum(cell_spk(j).t>times(i) & cell_spk(j).t<(times(i)+0.25))/0.25;
% end
% end
% 
% 
% figure
% %plot(squeeze(AllTrials.MUTcramp(10,1,1:138)), squeeze(spikerateCS(1,:)), '.');
% xlabel('CRamp')
% ylabel('spikerate during CS')
% %%
% 
% 
% R=corr(N1, AllTrials.MUTmean(10,round(linspace(1,200,40))));
% 
% figure
% plot(1:40, AllTrials.MUTmean(10,round(linspace(1,200,40))))
%% saving
Spikes(1).all_spikes = nan((size(cell_spk,2)-1)*139,10000);

for i=1:9 %%change this number according to the amount of cells, or the cells you want to save!!
    Spikes(i).per_cell=spikes(:,:,i);
    Spikes(1).all_spikes((139*i-138):i*139,1:10000)=spikes(:,:,i);
end

%%
save(strcat(path,'spikes.mat'),'Spikes')
print(psth, fullfile(path, 'psth-percell.png'), '-dpng')