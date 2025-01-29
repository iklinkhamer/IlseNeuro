function varargout = makePlots(trials,varargin)

ms2frm = @(ms) ms ./ 1e3 * 200;	% Convert ms to frames

if length(varargin) > 0
    isi = varargin{1};
    session = varargin{2};
    us = varargin{3};
    cs = varargin{4};
else
    isi = 200;
    session = 1; % Session s01
    us = 3;
    cs = 2;
end

%% Eyelid traces
hf1 = figure

hax=axes;

idx = find(trials.c_usnum==us & trials.c_csnum==cs & ismember(trials.session_of_day, session));

set(hax,'ColorOrder', jet(length(idx)), 'NextPlot', 'ReplaceChildren');
subplot(2,2,1)
plot(trials.tm(1,:), trials.eyelidpos(idx, :)')     % transpose matrix before plotting to deal with rare cases where m == n and Matlab chooses to plot columns instead of rows
hold on 
plot(trials.tm(1, :), mean(trials.eyelidpos(idx, :)), 'k', 'LineWidth', 2)
axis([trials.tm(1, 1) trials.tm(1, end) -0.1 1.1])
title('Conditioning: All data')
xlabel('Time from CS (s)')
ylabel('Eyelid pos (FEC)')

CS_US = trials.eyelidpos(trials.c_usdur~=0,:);
CSonly = trials.eyelidpos(trials.c_usdur==0,:);

subplot(2,2,2)
plot(trials.tm(1,1:41), CS_US(:,1:41), 'b', trials.tm(1,41:91),CS_US(:,41:91), 'r',trials.tm(1,91:199),CS_US(:,91:199),'b')
hold on
plot(trials.tm(1,:), mean(CS_US), 'k', 'LineWidth',2);
hold on
plot(trials.tm(1,41:96), -0.2*ones(56), 'c', trials.tm(1,91:96), -0.18*ones(6), 'g', 'LineWidth', 5);
ylim([min(min(trials.eyelidpos)) max(max(trials.eyelidpos))])
title('Conditioning: CS-US')
xlabel('Time from CS (s)')
ylabel('Eyelid pos (FEC)')


subplot(2,2,3)
plot(trials.tm(1,1:41), CSonly(:,1:41), 'b', trials.tm(1,41:91),CSonly(:,41:91), 'r',trials.tm(1,91:199),CSonly(:,91:199),'b')
hold on
plot(trials.tm(1,:), mean(CSonly), 'k', 'LineWidth',2);
hold on
plot(trials.tm(1,41:96), -0.2*ones(56), 'c', trials.tm(1,91:96), -0.18*ones(6), 'g', 'LineWidth', 5);
ylim([min(min(trials.eyelidpos)) max(max(trials.eyelidpos))])
title('Conditioning: CS only')
xlabel('Time from CS (s)')
ylabel('Eyelid pos')
%% CR amplitudes
pre = 1:ms2frm(200)
win = ms2frm(200):ms2frm(200 + isi)
cramp = mean(trials.eyelidpos(:, win), 2) - mean(trials.eyelidpos(:, pre), 2);
crabs = mean(trials.eyelidpos(:, win), 2);

subplot(2,2,4)

idx = find(trials.c_usnum==us & trials.c_csnum==cs & ismember(trials.session_of_day, session));
plot(trials.trialnum(idx), cramp(idx), '.')
hold on
plot(trials.trialnum(idx), crabs(idx), '+r')
plot([1 length(trials.trialnum)], [0.1 0.1], ':k')
axis([1 length(trials.trialnum) -0.1 1.1])
title('CS-US')
xlabel('Trials')
ylabel('CR size')

 if nargout > 0
    varargout{1} = hf1;
%     varargout{2} = hf2;
 end