%% IK 8-4-24
function plotEyelidTraces(mouseData, ID, fig, subcoords, day)
P = IkUtils.getParams();

lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(mouseData.onsets), 2)));
axes(fig);
set(gcf, 'DefaultFigureRenderer', 'painters');
subplot(subcoords(1), subcoords(2), subcoords(3))
hold on
if nargin < nargin('plotEyelidTraces')
    day = lastDay;
end

plot(mouseData.ts(1:41), squeeze(mouseData.traces(day,mouseData.type(:,1:P.nTrials)==1,1:41)),'color', [0.5 0.5 0.5])
plot(mouseData.ts(41:91), squeeze(mouseData.traces(day,mouseData.type(:,1:P.nTrials)==1,41:91)),'color', [0 0 1])
plot(mouseData.ts(91:95), squeeze(mouseData.traces(day,mouseData.type(:,1:P.nTrials)==1,91:95)),'color', [0 1 0])
plot(mouseData.ts(95:P.nTimeSteps), squeeze(mouseData.traces(day,mouseData.type(:,1:P.nTrials)==1,95:P.nTimeSteps)),'color', [0.5 0.5 0.5])
plot(mouseData.ts, nanmean(squeeze(mouseData.traces(day,mouseData.type(:,1:P.nTrials)==1,:))), 'k', 'LineWidth',2);

if length(ID) > 1
    Gene = ID(1);
    loc = ID(2);
    Type = ID(3);
    title(sprintf('All traces day %d: %s %s %s', day, Gene, loc, Type))
else
    title(sprintf('All traces day %d: %s', day, ID))
end

xlabel('Time from CS (ms)')
ylabel('fraction eyelid closure')
set(gca, "TickDir", "out")
box("off")
end