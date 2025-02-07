clear all

base_dir = 'D:\Nynke\neuroblink-data';

% mice = char('MI1915407','MI1915408', 'MI1930801', 'MI1960103');
% days = char('210510', '210511', '210512', '210513', '210514', '210517', '210518', '210519', '210520', '210521'); %Batch 1
%batch 2
mice = char('MI1944208'); %'MI1944203', 'MI1944205', 'MI1944206', 'MI1944208');
days = char('210517', '210518', '210519', '210520', '210521', '210524', '210525', '210526', '210527', '210531', '210601'); %Batch 1

%%
%days = char('210510', '210521'); %testing
%days = char('210510', '210511', '210512', '210513', '210514', '210517', '210518', '210519', '210520', '210521'); %Batch 1
%days = 1;%Batch 2
%sessions = ones(size(days,1), 1);
isi = 250;
us = 3;
cs = 1;
%CRabs = zeros(500);
%CRamp = zeros(500);
tic
for j = 1: size(mice,1)
    mouse = mice(j,:);
for i = 1:size(days,1)
day = days(i,:);
%session = sessions(i);
folder = fullfile(base_dir,  mouse, day);

trials(i) = processTrials(folder); %,'recalibrate');  % Recalibrate eyelid
folder2 = fullfile(base_dir,  mouse)
%data
save(fullfile(folder2, 'trialdata.mat'), 'trials');
% meanAll(i,:) = trials.meanAll;
% meanCS_US(i,:) = trials.meanCS_US;
% meanCSonly(i,:) = trials.meanCSonly;
% meanUSonly(i,:) = trials.meanUSonly;
% CRabs(i,1:length(trials.CRabs)) = trials.CRabs';
% CRamp(i,1:length(trials.CRamp)) = trials.CRamp';
end
time(j)=toc/60;
end
%% Plots

meanTraces = figure;

for i = 1:size(days,1)   
subplot(2,2,1)
hold on
plot(trials(i).tm(1,:), trials(i).meanAll)
axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
title('Conditioning: All data')
xlabel('Time from CS (ms)')
ylabel('mean eyelid pos')
legend(days)

subplot(2,2,2)
hold on
plot(trials(i).tm(1,:),trials(i).meanCS_US)
axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
title('Conditioning: CS-US')
xlabel('Time from CS (ms)')
ylabel('mean eyelid pos')

subplot(2,2,3)
hold on
plot(trials(i).tm(1,:), trials(i).meanCSonly)
axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
title('Conditioning: CS only')
xlabel('Time from CS (ms)')
ylabel('mean eyelid pos')

subplot(2,2,4)
hold on
plot(trials(i).tm(1,:), trials(i).meanUSonly)
axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
title('Conditioning: US only')
xlabel('Time from CS (ms)')
ylabel('mean eyelid pos')
end

alltraces = figure;
for i = 1:size(days,1)
% subplot(size(days,1),2,(2*i-1))
% plot(trials(i).tm(1,:), trials(i).eyelidpos');
% hold on
% plot(trials(i).tm(1,:), nanmean(trials(i).eyelidpos), 'k', 'LineWidth',2)   
% axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
% title( strcat(days(i,:), ': all traces - no processing'));
% xlabel('Time from CS (ms)')
% ylabel('eyelid pos')
subplot(size(days,1)/2,2,i)
plot(trials(i).tm(1,1:41), trials(i).tracesnorm(:,1:41), 'b', trials(i).tm(1,41:91),trials(i).tracesnorm(:,41:91), 'r',trials(i).tm(1,91:199),trials(i).tracesnorm(:,91:199),'b')
hold on
plot(trials(i).tm(1,:), nanmean(trials(i).tracesnorm), 'k', 'LineWidth',2);
axis([trials(i).tm(1, 1) trials(i).tm(1, end) min(min(trials(i).tracesnorm)) max(max(trials(i).tracesnorm))])
title( strcat(days(i,:), ': all traces - corrected'));
xlabel('Time from CS (ms)')
ylabel('eyelid pos')
end


allCRs = figure;

for i =1:1:size(days,1)
subplot(size(days,1),2,2*i-1)
plot(1:size(trials(i).CRamp), trials(i).CRamp, '.')
hold on
plot(1:size(trials(i).CRamp), 0.05*ones(size(trials(i).CRamp)))
title( strcat(days(i,:), ': all CRs'));
ylim([0 1])
subplot(size(days,1),2, 2*i)
plot(1:size(trials(i).CRamp(trials(i).CRamp>0.05)), trials(i).CRamp(trials(i).CRamp>0.05), '.')
hold on
plot(1:size(trials(i).CRamp(trials(i).CRamp>0.05)), 0.05*ones(size(trials(i).CRamp(trials(i).CRamp>0.05))))
title( strcat(days(i,:), ': CRs >5%'));
ylim([0 1])
end

for i =1:1:size(days,1)
ampCR(i) = trials(i).meanCR;
stdCR(i) = trials(i).stdCR;
ampCR5(i) = trials(i).meanCR5;
stdCR5(i) = trials(i).stdCR5;
end

meanCR = figure;
hold on
errorbar(ampCR,stdCR);
errorbar(ampCR5,stdCR5); 
xticks(1:size(days,1));
ylim([-0.2 1]);
a = get(gca,'XTickLabel');  
set(gca,'TickDir','out');
xlabel('Session');
ylabel('CR amplitude')
title('CS-US & CS-only trials')
box('off'); 
legend('CR amp', 'CR amp > 5% of UR mean')

%%
%check
figure
plot(trials(i).tm(1,:), trials(i).eyelidpos(trials(i).c_csdur==0,:));
%% saving
folder2 = fullfile(base_dir,  mouse)
%data
save(fullfile(folder2, 'trialdata.mat'), 'trials');

%figures
hgsave(meanTraces, fullfile(folder2, 'averageTraces.fig'));
hgsave(alltraces, fullfile(folder2, 'allTraces.fig'));
hgsave(allCRs, fullfile(folder2, 'allCRs.fig'));
hgsave(meanCR, fullfile(folder2, 'CRmean.fig'));

%pictures
print(meanTraces, fullfile(folder2, sprintf('%s_%s_meanTraces.png', mouse, day)), '-dpng')
print(alltraces, fullfile(folder2, sprintf('%s_%s_allTraces.png', mouse, day)), '-dpng')
print(allCRs, fullfile(folder2, sprintf('%s_%s_allCRs.png', mouse, day)), '-dpng')
print(meanCR, fullfile(folder2, sprintf('%s_%s_meanCRs.png', mouse, day)), '-dpng')


%% finding onset
% trace1 = trials(1).tracesnorm(180,:)
% plot(trials(1).tm(1,:), trace1)


%
% hgsave(hf1, fullfile(folder, 'meanTraces.fig'));
% print(hf1, fullfile(folder, sprintf('%s_%s_meanTraces.pdf', mouse, day)), '-dpdf')
% % 

% %% CR amplitudes
% 
% subplot(2,3,5)
% plot(1:256, CRamp(1,1:256), '.')
% % idx = find(trials.c_usnum==us & trials.c_csnum==cs & ismember(trials.session_of_day, session));
% % plot(trials.trialnum(idx), cramp(idx), '.')
% % hold on
% % plot(trials.trialnum(idx), crabs(idx), '+r')
% % plot([1 length(trials.trialnum)], [0.1 0.1], ':k')
% %axis([trials.trialnum(1) trials.trialnum(end) min(min(trials.eyelidpos)) max(max(trials.eyelidpos))])
% title('CRs')
% xlabel('Trials')
% ylabel('CR size')
% 
% 
% %%
% figure
% for i= [1 2 4 5 6 7]
%     subplot(2,4,i)
%     plot(trials.tm(1,:), trials.eyelidpos(trials.session_of_day == i,:))
%     ylim([-0.4 1.2])
% end