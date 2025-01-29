%load('trialdata.mat')
%% before processing
figure
subplot(2,2,1)
plot(trials.tm(1,:),trials.eyelidpos)
title('before processing')
%% put all traces to baseline
traces=bsxfun(@minus, trials.eyelidpos, median(trials.eyelidpos(:,1:40),2));

subplot(2,2,2)
plot(trials.tm(1,:), traces)
hold on
plot(trials.tm(1,:), nanmean(traces), 'k', 'LineWidth', 2)
title('baseline correction')
%% normalization
numSTD=3;
tracesnorm=nan(size(traces));


USmean = nanmedian(traces(trials.c_csdur ==0, 91:121), 2);
URmean=median(USmean);
tracesnorm = traces/URmean;
subplot(2,2,3)
plot(trials.tm(1,:), tracesnorm)
title('normalized with UR mean')
%% filter out invalid data

dat=tracesnorm(:,1:40);  %different time window, before US is presented (around 680 ms)
%dat =tracesnorm(idxses,t>0 & t<.4);
stds=nanstd(dat,[],2);%std of selection above, std of each trial
mus=median(stds); %median value of all the stds
iqrs=iqr(stds); %interquartile range (difference 75th - 25th percentiles of stds).
thresh=mus+numSTD*iqrs; %thresh = median of std + (3*iqrs-difference between interquartile range).
tracesnorm((stds>thresh),:)= nan ;
disp(['Excluded ' num2str(numel(find(stds>thresh))) ' Trials'])

subplot(2,2,4)
plot(trials.tm(1,:), tracesnorm)
hold on
plot(trials.tm(1,:), nanmean(tracesnorm), 'k', 'LineWidth', 2)
title('invalid data filtered out')

%% CR calculation

amplitudesCRs=nanmean(tracesnorm(:,75:90),2); %50 ms after CS
meanUR = nanmean(tracesnorm(:,91:121),2); %from 0 to 100 (100 ms aftr US)

CRs_5 = amplitudesCRs>(0.05*meanUR);
amplitudesCRs_5 = amplitudesCRs (amplitudesCRs>(0.05*meanUR));

amplitudeCRs_CSUS = grpstats(amplitudesCRs(trials.c_csdur==270 &trials.c_usdur==20), trials.session_of_day(trials.c_csdur==270 &trials.c_usdur==20), {'nanmean'})
std_amplitudeCRs_CSUS = grpstats(amplitudesCRs(trials.c_csdur==270 &trials.c_usdur==20), trials.session_of_day(trials.c_csdur==270 &trials.c_usdur==20), {'nanstd'})

amplitudeCRs_CSUS_5 = grpstats(amplitudesCRs(amplitudesCRs>(0.05*meanUR)&trials.c_csdur==270 &trials.c_usdur==20), trials.session_of_day(amplitudesCRs>(0.05*meanUR)&trials.c_csdur==270 &trials.c_usdur==20), {'nanmean'})
std_amplitudeCRs_CSUS_5 = grpstats(amplitudesCRs(amplitudesCRs>(0.05*meanUR)&trials.c_csdur==270 &trials.c_usdur==20), trials.session_of_day(amplitudesCRs>(0.05*meanUR)&trials.c_csdur==270 &trials.c_usdur==20), {'nanstd'})


figure
hold on
errorbar(amplitudeCRs_CSUS,std_amplitudeCRs_CSUS);
errorbar(amplitudeCRs_CSUS_5,std_amplitudeCRs_CSUS_5); 
% xlim([.5 10.5]);
% xticks([1:10]);
ylim([-0.5 1.5]);
a = get(gca,'XTickLabel');  
set(gca,'TickDir','out');
xlabel('Sessions');
ylabel('CR amplitude')
title('CS-US trials')
box('off'); 
legend('CR amp', 'CR amp > 5% of UR mean')
