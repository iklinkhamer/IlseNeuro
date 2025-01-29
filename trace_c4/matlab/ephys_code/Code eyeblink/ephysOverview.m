%% Nynke & I.K. 2-11-23
%WTpath= {'D:\Nynke\Ephys-data\MI19442-06F_2021-05-31_13-28-33_2\Record Node 101\'; 'D:\Nynke\Ephys-data\MI19442-06F_2021-06-01_12-40-32_3\Record Node 101\'};

%MUTpath = {'D:\Nynke\Ephys-data\MI19442-03F_2021-05-31_09-52-35_2\Record Node 101', 'D:\Nynke\Ephys-data\MI19442-03F_2021-06-01_10-15-12_3\Record Node 101', 'D:\Nynke\Ephys-data\MI19442-05F_2021-05-31_11-17-04_2\Record Node 101', 'D:\Nynke\Ephys-data\MI19442-05F_2021-06-01_11-26-56_3\Record Node 101', 'D:\Nynke\Ephys-data\MI19442-08M_2021-06-01_13-56-09_3\Record Node 101'};
function ephysOverview(mousenames, s)
arguments
    mousenames = mouseList(); 
    s = 1;
end
mousenamesWT = mousenames(1:length(mousenames)/2); % IK note: This list is absolutely not accurate, fix later 
mousenamesMT = mousenames(length(mousenames)/2+1:end);
n_mice = length(mousenames);
p = IkUtils.getParams();


data = struct;
for m = 1:n_mice
    mname = mousenames(m);
    path_extention = fullfile(mname,p.s(s));
    spikes_ = load(fullfile(p.pathSpikeSortingHome, path_extention, 'spikes.mat'));
    data(m).spikes = spikes_.spikes;
    trialdata_ = load(fullfile(p.pathSpikeSortingHome, path_extention, 'trialdata.mat'));
    data(m).trials = trialdata_.behavior_trial_data;

end
data.WTmask = zeros(1,n_mice);
data.WTmask(1:n_mice/2) = 1;
% 
% for i = 1:2
%     WTAllData(i)=load(fullfile(char(WTpath(i)), 'spikes.mat'));
% end
% for i = 1:5
%     MUTAllData(i)=load(fullfile(char(MUTpath(i)), 'spikes.mat'));
% end
%%

% load('D:\Nynke\neuroblink-data\MI1944206\trialdata.mat');

%%
% load('D:\Nynke\Ephys-data\MI19442-06F_2021-05-31_13-28-33_2\Record Node 101\data_sel.mat');
timestamps_ = load(fullfile(p.pathEphysDataDATA, path_extention,'timestamps.mat'));
timestamps = timestamps_.timestamps;
% load('D:\Nynke\Ephys-data\MI19442-06F_2021-05-31_13-28-33_2\Record Node 101\stimtimes.mat');
stimtimes_ = load(fullfile(p.pathEphysDataDATA, path_extention, 'stimtimes.mat'));
stimtimes = stimtimes_.stimtimes;
%%
% load('D:\Nynke\neuroblink-data\MI1944203\trialdata.mat');
% MUTmean = [trials(11).eyelidpos; trials(12).eyelidpos];
% load('D:\Nynke\neuroblink-data\MI1944205\trialdata.mat');
% MUTmean = [MUTmean; trials(11).eyelidpos; trials(12).eyelidpos];
% load('D:\Nynke\neuroblink-data\MI1944208\trialdata.mat');
% MUTmean = mean([MUTmean; trials(11).eyelidpos],1);

%% preprocessing
WTAllData(1).Spikes(4) = []; %use to delete empty fields before continuing!
MUTAllData(5).Spikes(2).all_spikes = MUTAllData(5).Spikes(1).all_spikes;
MUTAllData(5).Spikes(1) = [];
%% rasterplotting (turn off when not on MOAC!)
% kk=1;
% for i =1:2
%     for j=1:size(WTAllData(i).Spikes,2)
%         colour = [rand, rand, rand]; 
%         for k=1:138
%             hold on
%             plot(WTAllData(i).Spikes(j).per_cell(k,1:1000), (138*kk-k), '.', 'Color', colour);
%         end
%         kk=kk+1
%     end
% end
% xlim([-0.2 0.8])
% ylim([0 (kk-1)*138])
%
% kk=1;
% for i =1:5
%     for j=1:size(MUTAllData(i).Spikes,2)
%         colour = [rand, rand, rand]; 
%         for k=1:133
%             hold on
%             plot(MUTAllData(i).Spikes(j).per_cell(k,1:1000), (133*kk-k), '.', 'Color', colour);
%         end
%         kk=kk+1
%     end
% end
% xlim([-0.2 0.8])
% ylim([0 (kk-1)*133])


WTpsth = figure('Color', 'white','Position', [10 0 1000 1000]);

hold on
title('WT psth')
[N1,edges1]=histcounts([WTAllData(1).Spikes(1).all_spikes; WTAllData(2).Spikes(1).all_spikes],100);
a1 = area([-0.2 0.8], [max(N1) max(N1)]);
a1.FaceColor = [0.5,0.5,0.5];
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1 = area([0.25 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,1,0];
h = histogram('BinCounts', N1, 'BinEdges', edges1);
h.FaceColor = [0.9,0.9,0.9];
h.FaceAlpha = 0.8;
plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
plot(linspace(-0.2,0.8,200), WTmean*1200+100, 'LineWidth', 3, 'Color', [0,0,0])
legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
ylim([0 max(N1)])
xlim([-0.2 2])

MUTpsth = figure('Color', 'white','Position', [10 0 1000 1000]);
hold on
title('MUT psth')
[N1,edges1]=histcounts([MUTAllData(1).Spikes(1).all_spikes; MUTAllData(2).Spikes(1).all_spikes;MUTAllData(3).Spikes(1).all_spikes;MUTAllData(4).Spikes(1).all_spikes;MUTAllData(5).Spikes(1).all_spikes],100);
a1 = area([-0.2 0.8], [max(N1) max(N1)]);
a1.FaceColor = [0.5,0.5,0.5];
a1 = area([0 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,0,1];
a1 = area([0.25 0.27], [max(N1) max(N1)]);
a1.FaceColor = [0,1,0];
h = histogram('BinCounts', N1, 'BinEdges', edges1);
h.FaceColor = [0.9,0.9,0.9];
h.FaceAlpha = 0.8;
plot(linspace(-0.200,0.8,5), ones(5,1)*0, 'LineWidth', 10,'Color', [0.5,0.5,0.5])
plot(linspace(0,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,0,1])
plot(linspace(0.25,0.27,5), ones(5,1)*0, 'LineWidth', 10, 'Color', [0,1,0])
plot(linspace(-0.2,0.8,200), MUTmean*1000+50, 'LineWidth', 3, 'Color', [0,0,0])
legend('recording window', 'CS presented', 'US presented','psth', 'Location', 'southeast')
ylim([0 max(N1)])
xlim([-0.2 2])

