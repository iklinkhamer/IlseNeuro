%% IK 28-5-24
% Called by function: saveSummaryFiguresMouseWrapperIK
function plotSummaryFiguresMUTvsWT_IK(mice, gene, loc, outputRoot, loadPrevData)
arguments
    mice
    gene
    loc
    outputRoot
    loadPrevData
end

if ~loadPrevData || ~exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")), 'file')
    % Initialize empty arrays for WT
    allStatsWT = struct;
    allStatsWT.cs = [];
    allStatsWT.us = [];

    % Populate the structures for WT
    for m = 1:length(mice{1})
        mouseCodesGroup = mice{1};
        allStatsWTtemp = histStatsMouseWrapper( ...
            mouseCodesGroup(m), ...
            'event', "cs_facilitation", ...
            'modulating', 1 ...
            );
        if ~isempty(allStatsWTtemp)
            allStatsWT.cs = [allStatsWT.cs allStatsWTtemp.cs];  % Concatenate 'cs' field
            allStatsWT.us = [allStatsWT.us allStatsWTtemp.us];  % Concatenate 'us' field
        end
    end

    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")), "allStatsWT");
else
    load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "WT")));
end

if ~loadPrevData || ~exist(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")), 'file')   
    allStatsMUT = struct;   % Initialize empty arrays for MUT
    allStatsMUT.cs = [];
    allStatsMUT.us = [];
    
    for m = 1:length(mice{2})% Populate the structures for MUT
        mouseCodesGroup = mice{2};
        allStatsMUTtemp = histStatsMouseWrapper( ...
            mouseCodesGroup(m), ...
            'event', "cs_facilitation", ...
            'modulating', 1 ...
            );
        if ~isempty(allStatsMUTtemp)
            allStatsMUT.cs = [allStatsMUT.cs allStatsMUTtemp.cs];  % Concatenate 'cs' field
            allStatsMUT.us = [allStatsMUT.us allStatsMUTtemp.us];  % Concatenate 'us' field
        end
    end

    save(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")), "allStatsMUT");
else
    load(fullfile("/home/i.klinkhamer/Documents/Data/other/",sprintf("allStatsBoxplots_%s%s%s.mat", gene, loc, "MUT")));
end

%% CS

%% covariance
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "cov");
makeBoxPlots(boxplot_data, boxplot_data_cell, "CV", gene, loc, outputRoot)

%% local covariance
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "meanLcov");
makeBoxPlots(boxplot_data, boxplot_data_cell, "LCV", gene, loc, outputRoot)

%% ISI
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "meanISIroi");
makeBoxPlots(boxplot_data, boxplot_data_cell, "ISI", gene, loc, outputRoot)

%% rate
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "rate");
makeBoxPlots(boxplot_data, boxplot_data_cell, "freq", gene, loc, outputRoot)

%% maxAmpTime
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "maxAmpTime");
makeBoxPlots(boxplot_data, boxplot_data_cell, "PT", gene, loc, outputRoot)

%% SD
[boxplot_data, boxplot_data_cell] = getBoxplotData([allStatsWT, allStatsMUT], "sd");
makeBoxPlots(boxplot_data, boxplot_data_cell, "SD", gene, loc, outputRoot)

%% Overall stats
figurePosition = [100, 100, 1800, 500]; figure('Position', figurePosition);
subplot(1,3,1)
hold on 
boxplot_data = [ ...
    [allStatsMUT.us.rateFull] ...
    [allStatsWT.us.rateFull] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.rateFull]; ...
    [allStatsWT.us.rateFull]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
positions = 1:length(boxplot_labels);
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
x_max = max(boxplot_data); % Maximum x-value of the data
x_offset = x_max*0.1; % Offset above the max value for the asterisks
for i = 1:(length(boxplot_labels)/2)
    idx1 = 2*i-1;
    idx2 = 2*i;
    if ~isempty(boxplot_data_cell{idx1}) && ~isempty(boxplot_data_cell{idx2})
        disp("Rate full window")
        [~, p, ~, stats] = ttest2(boxplot_data_cell{idx1},boxplot_data_cell{idx2})         
        if p < 0.001
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '***', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.01
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '**', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.05
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '*', 'FontSize', 10, 'VerticalAlignment', 'middle');
        else
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), 'n.s.', 'FontSize', 10, 'VerticalAlignment', 'middle');
        end
    end
end
ylabel('Type')
xlabel('frequency (Hz)')
xlim([0 250])
% xlim([min(boxplot_data), x_max + x_offset])
set(gca,'tickdir','out','box','off');
title(sprintf('Rate full trial window spk times MUT vs WT'))
hold off

subplot(1,3,2)
hold on 
boxplot_data = [ ...
    [allStatsMUT.us.baserate] ...
    [allStatsWT.us.baserate] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.baserate]; ...
    [allStatsWT.us.baserate]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
positions = 1:length(boxplot_labels);
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
x_max = max(boxplot_data); % Maximum x-value of the data
x_offset = x_max*0.1; % Offset above the max value for the asterisks
for i = 1:(length(boxplot_labels)/2)
    idx1 = 2*i-1;
    idx2 = 2*i;
    if ~isempty(boxplot_data_cell{idx1}) && ~isempty(boxplot_data_cell{idx2})
        disp("Baseline full window")
        [~, p, ~, stats] = ttest2(boxplot_data_cell{idx1},boxplot_data_cell{idx2})
        if p < 0.001
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '***', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.01
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '**', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.05
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '*', 'FontSize', 10, 'VerticalAlignment', 'middle');
        else
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), 'n.s.', 'FontSize', 10, 'VerticalAlignment', 'middle');
        end
    end
end
ylabel('Type')
xlabel('frequency (Hz)')
xlim([0 250])
% xlim([min(boxplot_data), x_max + x_offset])
set(gca,'tickdir','out','box','off');
title(sprintf('Baseline spk times MUT vs WT'))
hold off

subplot(1,3,3); hold on
boxplot_data = [ ...
    [allStatsMUT.us.covFull] ...
    [allStatsWT.us.covFull] ...
    ];
boxplot_data_cell = { ...
    [allStatsMUT.us.covFull]; ...
    [allStatsWT.us.covFull]; ...
    };
grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
boxplot_labels = {'MUT','WT'};
positions = 1:length(boxplot_labels);
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [[0.7 0 0]; [0 0 0.7];[0.7 0 0];[0 0 0.7]];
boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
x_max = max(boxplot_data); % Maximum x-value of the data
x_offset = x_max*0.1; % Offset above the max value for the asterisks
for i = 1:(length(boxplot_labels)/2)
    idx1 = 2*i-1;
    idx2 = 2*i;
    if ~isempty(boxplot_data_cell{idx1}) && ~isempty(boxplot_data_cell{idx2})
        disp("Covariance full window")
        [~, p, ~, stats] = ttest2(boxplot_data_cell{idx1},boxplot_data_cell{idx2})
        if p < 0.001
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '***', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.01
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '**', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p < 0.05
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '*', 'FontSize', 10, 'VerticalAlignment', 'middle');
        else
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), 'n.s.', 'FontSize', 10, 'VerticalAlignment', 'middle');
        end
    end
end
ylabel('Type')
xlabel('CV')
xlim([0 1]/10)

% xlim([min(boxplot_data), x_max + x_offset])
set(gca,'tickdir','out','box','off');
title(sprintf('Covariance spk times MUT vs WT'))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotall_rate_and_cov_MUTvsWT_%s%s.png', gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotall_rate_and_cov_MUTvsWT_%s%s.eps', gene, loc)))

end

function makeBoxPlots(boxplot_data, boxplot_data_cell, stat_name, gene, loc, outputRoot)
figurePosition = [100, 100, 1200, 500]; figure('Position', figurePosition);
hold on

boxplot_labels = {'nonMUT','nonWT','facMUT','facWT', 'supMUT', 'supWT'};
positions = 1:length(boxplot_labels);

subplot(1, 2, 1); hold on

boxplot_data_cs = boxplot_data{1};
boxplot_data_cell_cs = boxplot_data_cell(:,1);

grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell_cs{i}),1)},(1:numel(boxplot_data_cell_cs))'));
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [];
for l = 1:length(boxplot_labels_new)
    label = boxplot_labels_new{l};
    if contains(label,"MUT")
        colors = [colors; [0.7 0 0]];
    else
        colors = [colors; [0 0 0.7]];
    end
end
boxplot(boxplot_data_cs, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
x_max = max(boxplot_data_cs); % Maximum x-value of the data
x_offset = x_max*0.1; % Offset above the max value for the asterisks
for i = 1:(length(boxplot_labels)/2)
    idx1 = 2*i-1;
    idx2 = 2*i;
    if ~isempty(boxplot_data_cell_cs{idx1}) && ~isempty(boxplot_data_cell_cs{idx2})
        disp("CS")
        [~, p, ~, stats] = ttest2(boxplot_data_cell_cs{idx1},boxplot_data_cell_cs{idx2})
        if p <= 0.001
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '***', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p <= 0.01
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '**', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p <= 0.05
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '*', 'FontSize', 10, 'VerticalAlignment', 'middle');
        else
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), 'n.s.', 'FontSize', 10, 'VerticalAlignment', 'middle');
        end
    end
end
ylabel('Type')
xlabel(stat_name)
% xlim([min(boxplot_data_cs), x_max + x_offset])
if stat_name == "CV"
    xlim([0 5]/1000)
elseif stat_name == "LCV"
    xlim([1.5 2.5]/1000000)
elseif stat_name == "ISI"
    xlim([0 0.05])
elseif stat_name == "freq"
    xlim([0 250])
elseif stat_name == "PT"
    xlim([0 0.4])
elseif stat_name == "SD"
    xlim([25 75]/1000)
end
set(gca,'tickdir','out','box','off');
title(sprintf('%s spk times CS MUT vs WT', stat_name))
hold off

subplot(1,2,2); hold on

boxplot_data_us = boxplot_data{2};
boxplot_data_cell_us = boxplot_data_cell(:,2);

grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell_us{i}),1)},(1:numel(boxplot_data_cell_us))'));
boxplot_labels_new = {boxplot_labels{unique(grp)}};
colors = [];
for l = 1:length(boxplot_labels_new)
    label = boxplot_labels_new{l};
    if contains(label,"MUT")
        colors = [colors; [0.7 0 0]];
    else
        colors = [colors; [0 0 0.7]];
    end
end
boxplot(boxplot_data_us, grp, 'labels', boxplot_labels_new, 'orientation','horizontal','colors',colors)
x_max = max(boxplot_data_us); % Maximum x-value of the data
x_offset = x_max*0.1; % Offset above the max value for the asterisks
for i = 1:(length(boxplot_labels)/2)
    idx1 = 2*i-1;
    idx2 = 2*i;
    if ~isempty(boxplot_data_cell_us{idx1}) && ~isempty(boxplot_data_cell_us{idx2})
        disp("US")
        [h, p, ~, stats] = ttest2(boxplot_data_cell_us{idx1},boxplot_data_cell_us{idx2})
        if p <= 0.001
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '***', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p <= 0.01
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '**', 'FontSize', 10, 'VerticalAlignment', 'middle');
        elseif p <= 0.05
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), '*', 'FontSize', 10, 'VerticalAlignment', 'middle');
        else
            plot([x_max+x_offset, x_max+x_offset], [positions(idx1), positions(idx2)], 'k-', 'LineWidth', 1);
            text(x_max+x_offset,mean([positions(idx1), positions(idx2)]), 'n.s.', 'FontSize', 10, 'VerticalAlignment', 'middle');
        end
    end
end
ylabel('Type')
xlabel(stat_name)
if stat_name == "CV"
    xlim([0 5]/1000)
elseif stat_name == "LCV"
    xlim([1.5 2.5]/1000000)
elseif stat_name == "ISI"
    xlim([0 0.05])
elseif stat_name == "freq"
    xlim([0 250])
elseif stat_name == "PT"
    xlim([0 0.4])
elseif stat_name == "SD"
    xlim([25 75]/1000)
end
% xlim([min(boxplot_data_us), x_max + x_offset])
set(gca,'tickdir','out','box','off');
title(sprintf('%s spk times US MUT vs WT', stat_name))
hold off

saveas(gcf, fullfile(outputRoot, sprintf('boxplotall%s_MUTvsWT_%s%s.png', stat_name, gene, loc)), 'png')
print(gcf, '-depsc', '-painters', fullfile(outputRoot, sprintf('boxplotall%s_MUTvsWT_%s%s.eps', stat_name, gene, loc)))
end
