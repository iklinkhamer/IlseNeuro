%% I.K. 1-6-24
function plotOffset(mouseData, ID, fig, subcoords, day)
onsets_firstday_ = mouseData.onsets(1, :);% Extract onset times for WT and MUT
onsets_lastday_ = mouseData.onsets(day, :);

onsets_firstday = onsets_firstday_(~isnan(onsets_firstday_));% Remove NaN values
onsets_lastday = onsets_lastday_(~isnan(onsets_lastday_));

[h, p] = kstest2(onsets_firstday, onsets_lastday);% Perform KS test
disp(['KS test result: h = ', num2str(h), ', p = ', num2str(p)])% Display KS test result

% Perform Mann-Whitney U test
[p_mw, h_mw] = ranksum(onsets_firstday, onsets_lastday);
disp(['Mann-Whitney U test result: h = ', num2str(h_mw), ', p = ', num2str(p_mw)])

% Perform Wilcoxon Signed-Rank Test
onsets_firstday = onsets_firstday_(~isnan(onsets_firstday_) & ~isnan(onsets_lastday_));
onsets_lastday = onsets_lastday_(~isnan(onsets_firstday_) & ~isnan(onsets_lastday_));
[p_ws, h_ws] = signrank(onsets_firstday, onsets_lastday);
disp(['Wilcoxon Signed-Rank Test: h = ', num2str(h_ws), ', p = ', num2str(p_ws)]);


num_trials = size(mouseData.onsets, 2); % Total number of wildtype trials
lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(mouseData.onsets), 2)));
if nargin < 5
    day = lastDay;
end
axes(fig);
set(gcf, 'DefaultFigureRenderer', 'painters');
subplot(subcoords(1), subcoords(2), subcoords(3))
hold on

binSize = 25; 
binEdges = 0:binSize:250;
barOffset = binSize/2;

% CR onset time for day 1
wt_counts = histcounts(mouseData.onsets(1,:), binEdges);
wt_hist_norm = wt_counts / num_trials;
bar(binEdges(1:end-1) + barOffset, wt_hist_norm, 'FaceColor', 'b', 'EdgeColor', 'none', 'BarWidth', 1, 'FaceAlpha', 0.5);

% CR onset time for last day
wt_counts = histcounts(mouseData.onsets(day,:), binEdges);
wt_hist_norm = wt_counts / num_trials;
bar(binEdges(1:end-1)+barOffset, wt_hist_norm, 'FaceColor', 'r', 'EdgeColor', 'none', 'BarWidth', 1, 'FaceAlpha', 0.5);

xlim([0 250])
ylim([0 0.3])
xticks(binEdges);
if length(ID) > 1
    Gene = ID(1);
    loc = ID(2);
    Type = ID(3);
    title(sprintf('%s %s %s CR Onset Latency', Gene, loc, Type))
else
    title(sprintf('%s CR Onset Latency', ID))
end

if subcoords(3) == 1 || subcoords(3) == 3
    legend('Day1', sprintf("Day%d", day))
end
xlabel('Time (ms)')
ylabel('Normalized Counts per Trial')
set(gca,'TickDir','out');
box('off');
end