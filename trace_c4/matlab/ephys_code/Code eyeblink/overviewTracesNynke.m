clear all
close all
%% preallocation

base_dir = 'D:\Nynke\neuroblink-data';

MUTmice = char('MI1915407','MI1915408', 'MI1930801','MI1944203', 'MI1944205', 'MI1944208' );
WTmice = char('MI1944206','MI1960103' );

type = ones(1,240);
WTtraces=nan(10,size(WTmice,1)*240,200);
% WTtype=ones(10,size(WTmice,1)*240);
WTcramp=nan(10,size(WTmice,1),220);
WTcramp5=nan(10,size(WTmice,1),220);
WTmean=nan(10,200);
WTonsets=nan(10,size(WTmice,1)*240);
MUTtraces=nan(10,size(MUTmice,1)*240,200);
% MUTtype=ones(10,size(MUTmice,1)*240);
MUTcramp=nan(10,size(MUTmice,1),220);
MUTcramp5=nan(10,size(MUTmice,1),220);
MUTmean=nan(10,200);
MUTonsets=nan(10,size(MUTmice,1)*240);
%% calculations

for i=1:size(WTmice,1)
    mouse = WTmice(i,:);
    folder = fullfile(base_dir,  mouse);
    load(strcat(folder, '\trialdata.mat'));
    for j=1:10
        WTtraces(j,(i*240-239):i*240,:) = trials(j).tracesnorm(1:240,:);
        WTcramp(j,i,1:220) = trials(j).CRamp(1:220)';
        WTcramp5(j,i,1:length(trials(j).CRamp5)) = trials(j).CRamp5';
        WTCRperc(j,i) = sum(~isnan(WTcramp5(j,i,:)))/220*100;
    end
end


for i=1:size(MUTmice,1)
    mouse = MUTmice(i,:);
    folder = fullfile(base_dir,  mouse);
    load(strcat(folder, '\trialdata.mat'));
    for j=1:10
        MUTtraces(j,(i*240-239):i*240,:) = trials(j).tracesnorm(1:240,:);
        MUTcramp(j,i,1:220) = trials(j).CRamp(1:220)';
        MUTcramp5(j,i,1:length(trials(j).CRamp5)) = trials(j).CRamp5';
        MUTCRperc(j,i) = sum(~isnan(MUTcramp5(j,i,:)))/220*100;
    end
end


type(trials(1).c_csdur == 0) = 0;
type(trials(1).c_usdur == 0) = 2;
WTtype = [ type type];
MUTtype = [ type type type type type type];

for j =1:10
    WTcrampmean(j)=nanmean(WTcramp(j,:,:), 'all');
    WTcramp5mean(j)=nanmean(WTcramp5(j,:,:), 'all');
    WTcrampstd(j) = nanstd(WTcramp(j,:,:),0, 'all');
    WTcramp5std(j) = nanstd(WTcramp5(j,:,:),0,'all');
    WTmean(j,:) = nanmean(WTtraces(j,WTtype == 1,:),2);
    
    MUTcrampmean(j)=nanmean(MUTcramp(j,:,:), 'all');
    MUTcramp5mean(j)=nanmean(MUTcramp5(j,:,:), 'all');
    MUTcrampstd(j) = nanstd(MUTcramp(j,:,:),0, 'all');
    MUTcramp5std(j) = nanstd(MUTcramp5(j,:,:),0,'all');
    MUTmean(j,:) = nanmean(MUTtraces(j,MUTtype == 1,:),2);
end

WTCRperdstd = std(WTCRperc,1,2);
MUTCRperdstd = std(MUTCRperc,1,2);
ts=trials(1).tm(1,:);

%CR onset
for i = 1:10
    for j=1:(length(WTtraces))
        if ~isnan(WTtraces(i,j,1))
            baseline = nanmean(WTtraces(i,j,1:40));
            stdbase = nanstd(WTtraces(i,j,1:40));
            
            if find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1)
                WTonsets(i,j) = ts(find(WTtraces(i,j,:)>(baseline+5*stdbase) , 1));
            end
        end
    end
end

WTonsets = WTonsets(:,WTtype ==1 |WTtype ==2);
%WTonsets(WTcramp <0.05) = NaN;

for i = 1:10
    for j=1:(length(MUTtraces))
        if ~isnan(MUTtraces(i,j,1))
            baseline = nanmean(MUTtraces(i,j,1:40));
            stdbase = nanstd(MUTtraces(i,j,1:40));
            
            if find(MUTtraces(i,j,:)>(baseline+10*stdbase) , 1)>0
                MUTonsets(i,j) = ts(find(MUTtraces(i,j,:)>(baseline+10*stdbase) , 1));
            end
        end
    end
end

MUTonsets = MUTonsets(:,MUTtype ==1 |MUTtype ==2);
%MUTonsets(MUTcramp <0.05) = NaN;
% %% plotting
% meanTraces = figure('Color', 'white', 'Position', [0 0 1500 2000]);
% subplot(1,2,1)
% plot(ts, WTmean)
% legend('day 1', 'day 2','day 3','day 4','day 5','day 6','day 7','day 8','day 9','day 10')
% axis([ts(1) ts(end) min(min(WTmean)) max(max(WTmean))])
% title('Mean eyelid traces - WT');
% xlabel('Time from CS (ms)')
% ylabel('fraction eyelid closure')
% subplot(1,2,2)
% plot(ts, MUTmean)
% legend('day 1', 'day 2','day 3','day 4','day 5','day 6','day 7','day 8','day 9','day 10')
% axis([ts(1) ts(end) min(min(MUTmean)) max(max(MUTmean))])
% title('Mean eyelid traces - MUT');
% xlabel('Time from CS (ms)')
% ylabel('fraction eyelid closure')
% 
% CRamp = figure('Color', 'white');
% hold on
% errorbar(WTcramp5mean,WTcramp5std, 'b');
% errorbar(MUTcramp5mean,MUTcramp5std, 'r'); 
% xticks(1:10);
% ylim([-0.2 1]);
% a = get(gca,'XTickLabel');  
% set(gca,'TickDir','out');
% xlabel('Session');
% ylabel('CR amplitude')
% title('CS-US & CS-only trials - CR>5%')
% box('off'); 
% legend('WT', 'MUT')
% 
% CRamp5 = figure('Color', 'white');
% hold on
% errorbar(WTcrampmean,WTcrampstd, 'b');
% errorbar(MUTcrampmean,MUTcrampstd, 'r'); 
% xticks(1:10);
% ylim([-0.2 1]);
% a = get(gca,'XTickLabel');  
% set(gca,'TickDir','out');
% xlabel('Session');
% ylabel('CR amplitude')
% title('CS-US & CS-only trials ')
% box('off'); 
% legend('WT', 'MUT')

%% Final figures
eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

subplot(2,2,1)
hold on
plot(ts(1:41), squeeze(WTtraces(10,WTtype(1,1:240)==1,1:41)),'color', [0.5 0.5 0.5])
plot(ts(41:91), squeeze(WTtraces(10,WTtype(1,1:240)==1,41:91)),'color', [0 0 1])
plot(ts(91:95), squeeze(WTtraces(10,WTtype(1,1:240)==1,91:95)),'color', [0 1 0])
plot(ts(95:200), squeeze(WTtraces(10,WTtype(1,1:240)==1,95:200)),'color', [0.5 0.5 0.5])
plot(ts, nanmean(squeeze(WTtraces(10,WTtype(1,1:240)==1,:))), 'k', 'LineWidth',2);
title('All traces day 10: WT')
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,2)
hold on
MUTtracesPlot = MUTtraces(:,481:720,:);
plot(ts(1:41), squeeze(MUTtracesPlot(10,MUTtype(1,1:240)==1,1:41)),'color', [0.5 0.5 0.5])
plot(ts(41:91), squeeze(MUTtracesPlot(10,MUTtype(1,1:240)==1,41:91)),'color', [0 0 1])
plot(ts(91:95), squeeze(MUTtracesPlot(10,MUTtype(1,1:240)==1,91:95)),'color', [0 1 0])
plot(ts(95:200), squeeze(MUTtracesPlot(10,MUTtype(1,1:240)==1,95:200)),'color', [0.5 0.5 0.5])
plot(ts, nanmean(squeeze(MUTtracesPlot(10,MUTtype(1,1:240)==1,:))), 'k', 'LineWidth',2);
title('All traces day 10: MUT')
xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')

subplot(2,2,3)
WTmeanWF=WTmean;
WTmeanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
WTmeanWF(:,200)=NaN;
colormap(winter)
for i=1:10
color(i,1:200)=i*(length(WTmeanWF)); %give each line a different color (instead of gradient for heigth profile)
end
h1=waterfall(ts,1:10,WTmeanWF,color)
h1.LineWidth=2;
h1.FaceAlpha=0; %transparent faces
view([-0.5 -1.5 1]) %determines view angle
zlim=([0 1.2]);
yticks(1:10);
xlim([-200 800]);
title('Average traces day 1-10: WT')
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')

subplot(2,2,4)
MUTmeanWF=MUTmean;
MUTmeanWF(:,1)=NaN;
MUTmeanWF(:,200)=NaN;
h2=waterfall(ts,1:10,MUTmeanWF,color)
h2.LineWidth=2;
h2.FaceAlpha=0;
view([-0.5 -1.5 1])
zlim=([0 1.2]);
yticks(1:10);
title('Average traces day 1-10: MUT')
xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel('training day')


CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

subplot(2,2,1)
hold on
errorbar(WTcrampmean,WTcrampstd, 'b');
errorbar(MUTcrampmean,MUTcrampstd, 'r'); 
xticks(1:10);
ylim([-0.2 1]);
a = get(gca,'XTickLabel');  
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR amplitude')
title('CS-US & CS-only trials')
box('off'); 
legend('WT', 'MUT')

subplot(2,2,2)
hold on
errorbar(mean(WTCRperc,2),std(WTCRperc,1,2), 'b');
errorbar(mean(MUTCRperc,2),std(MUTCRperc,1,2), 'r'); 
xticks(1:10);
ylim([0 100]);
a = get(gca,'XTickLabel');  
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR percentage')
title('CS-US & CS-only trials')
box('off'); 

subplot(2,2,3)
%CR onset time
histogram(WTonsets(1,:))
hold on
histogram(WTonsets(9,:))
xlim([0 350])
title('WT onset latency')
legend('day 1', 'day10')

subplot(2,2,4)
%heatmap/density plot
histogram(MUTonsets(1,:))
hold on
histogram(MUTonsets(9,:)) % right now day 9 because better results...
xlim([0 350])
title('MUT onset latency')


% %% saving
% AllTrials.WTtraces=WTtraces;
% AllTrials.MUTtraces=MUTtraces;
% AllTrials.WTmean=WTmean;
% AllTrials.MUTmean=MUTmean;
% AllTrials.WTcramp=WTcramp;
% AllTrials.WTcrampstd=WTcrampstd;
% AllTrials.MUTcramp=MUTcramp;
% AllTrials.MUTcrampstd=MUTcrampstd;
% AllTrials.WTcramp5=WTcramp5;
% AllTrials.WTcramp5std=WTcramp5std;
% AllTrials.MUTcramp5=MUTcramp5;
% AllTrials.MUTcramp5std=MUTcramp5std;
% 
% 
% save(fullfile(base_dir, 'AllTrials.mat'), 'AllTrials');
% % 
% %figures
% hgsave(meanTraces, fullfile(base_dir, 'meanTraces.fig'));
% hgsave(eyeblinkTracesOverview, fullfile(base_dir, 'eyeblinkTracesOverview.fig'));
% hgsave(CRamp, fullfile(base_dir, 'CRamp.fig'));
% hgsave(CRamp5, fullfile(base_dir, 'CRamp5.fig'));
% 
% %pictures
% print(meanTraces, fullfile(base_dir, 'meanTraces.png'), '-dpng')
% print(eyeblinkTracesOverview, fullfile(base_dir, 'eyeblinkTracesOverview.png'), '-dpng')
% print(CRamp, fullfile(base_dir, 'CRamp.png'), '-dpng')
% print(CRamp5, fullfile(base_dir, 'CRamp5.png'), '-dpng')

print(eyeblinkTracesOverview, fullfile(base_dir, 'eyeblinkTracesOverview.png'), '-dpng')
print(CROverview, fullfile(base_dir, 'CROverview.png'), '-dpng')