baseline = mean(mean(trials.eyelidpos(trials.c_csdur==0,1:40)));
fullclosure = mean(max(trials.eyelidpos(trials.c_csdur==0,96:200),[],2));
trials.traces = (trials.eyelidpos - baseline) / (fullclosure - baseline);
trials.meanAll = mean(trials.traces);
trials.meanCS_US = mean(trials.traces(trials.c_usdur~=0,:));
trials.meanCSonly = mean(trials.traces(trials.c_usdur==0,:));
trials.meanUSonly = mean(trials.traces(trials.c_csdur~=0,:));

%%
figure
subplot(2,2,1)
plot(trials.tm(1,:), trials.eyelidpos)
% hold on
% plot(trials.tm(1,:), trials.meanAll, 'k', 'LineWidth', 2)
axis([trials.tm(1, 1) trials.tm(1, end) min(min(trials.eyelidpos)) max(max(trials.eyelidpos))])
title('All data')
xlabel('Time from CS (s)')
ylabel('eyelid pos')
hold on
plot(trials.tm(1,41:96), -0.2*ones(56), 'c', trials.tm(1,91:96), -0.18*ones(6), 'g', 'LineWidth', 5);

subplot(2,2,2)
plot(trials.tm(1,:), trials.traces(trials.c_usdur~=0,:))
hold on
plot(trials.tm(1,:), trials.meanCS_US, 'k', 'LineWidth', 2)
title('CS-US recalibrated')
hold on
plot(trials.tm(1,41:96), -0.2*ones(56), 'c', trials.tm(1,91:96), -0.18*ones(6), 'g', 'LineWidth', 5);

subplot(2,2,3)
plot(trials.tm(1,:), trials.traces(trials.c_usdur==0,:))
hold on
plot(trials.tm(1,:), trials.meanCSonly, 'k', 'LineWidth', 2)
title('CS-only recalibrated')
hold on
plot(trials.tm(1,41:96), -0.2*ones(56), 'c', trials.tm(1,91:96), -0.18*ones(6), 'g', 'LineWidth', 5);

subplot(2,2,4)
plot(trials.trialnum, CRamp(15,1:240), '.')
hold on
%plot(trials.trialnum, CRabs(15,1:240), '+r')
plot([1 length(trials.trialnum)], [0.1 0.1], ':k')
axis([trials.trialnum(1) trials.trialnum(end) min(min(trials.eyelidpos)) max(max(trials.eyelidpos))])
%legend('CR amp', 'CR abs')
title('CR amplitude')


