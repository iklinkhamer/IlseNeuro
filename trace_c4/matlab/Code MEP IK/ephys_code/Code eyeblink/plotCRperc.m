%% IK 8-4-24
function plotCRperc(Data, ID, fig, subcoords, kwargs)
arguments
    Data = [];
    ID = [];
    fig = [];
    subcoords = [1 1 1];
    kwargs.showBackGroundLines = true;
    kwargs.day = min(IkUtils.getParams().lastDay, sum(any(~isnan(Data.onsets), 2)));
end
showBackGroundLines = kwargs.showBackGroundLines;
day = kwargs.day;
if length(Data) < 2
%     lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(Data.onsets), 2)));
%     if nargin < 6
%         day = lastDay;
%     end
    axes(fig); hold on
    % set(gcf, 'DefaultFigureRenderer', 'painters');
    subplot(subcoords(1), subcoords(2), subcoords(3))
    hold on
    plot(1:day,Data.CRperc(1:day), 'b','LineWidth',1.2);
    xticks(1:day);
    xlim([1 day]);
    ylim([0 100]);
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR percentage')
    title('CR percentage during training')
    box('off');

else
    WT = Data(1);
    MUT = Data(2);
%     lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(WT.onsets), 2)));
%     if nargin < 6
%         day = lastDay;
%     end

    axes(fig);
    % set(gcf, 'DefaultFigureRenderer', 'painters');
    subplot(subcoords(1), subcoords(2), subcoords(3))
    hold on
    if showBackGroundLines
        plot(WT.CRperc(:,:,1), 'color', [0.68, 0.85, 0.9, 0.5]);
        plot(MUT.CRperc(:,:,1), 'color', [1, 0.8,0.8, 0.5]);
    end
    errorbar(nanmean(WT.CRperc,2),WT.CRpercSEM, 'b','LineWidth',1.2);
    errorbar(nanmean(MUT.CRperc,2),MUT.CRpercSEM, 'r','LineWidth',1.2);
    xticks(1:day);
    xlim([1, day]);
    ylim([0 100]);
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR percentage')
    if length(ID) > 1
        Gene = ID(1);
        loc = ID(2);
        title(sprintf('CR percentage %s %s', Gene, loc))
    else
        title(sprintf('CR percentage %s', ID))
    end
    
    box('off');

    repmAnova(WT,MUT, day, Gene, loc)
end
end

function repmAnova(WT,MUT, day, Gene, loc)
P = IkUtils.getParams();
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
Sessions = table((1:day)', 'VariableNames', {'Sessions'});

WT_CRperc = WT.CRperc(1:day,:);
MUT_CRperc = MUT.CRperc(1:day,:);
[~,colWT] = find(isnan(WT_CRperc));
WT_CRperc(:, colWT) = [];
[~,colMUT] = find(isnan(MUT_CRperc));
MUT_CRperc(:,colMUT) = [];

variable_names = [{'Group'}, strcat('s', string(1:day))];
variable_names_cell = cellfun(@(x) {x}, variable_names);

[sessionsWT, miceWT, trialsWT] = size(WT_CRperc);
reshaped_data_WT_new = reshape(WT_CRperc, [sessionsWT, miceWT*trialsWT]);
groupWT = repmat({'WT'}, size(reshaped_data_WT_new, 2),1);

[sessionsMUT, miceMUT, trialsMUT] = size(MUT_CRperc);
reshaped_data_MUT = reshape(MUT_CRperc, [sessionsMUT, miceMUT*trialsMUT]);
groupMUT = repmat({'MUT'}, size(reshaped_data_MUT, 2),1);

concatenated_array = [reshaped_data_WT_new, reshaped_data_MUT]';
group = [groupWT;groupMUT];

combined_cell = cell(size(concatenated_array,1), day+1);
combined_cell(:, 1:size(group, 2)) = group; % Fill the cell array with the data
combined_cell(:, size(group, 2)+1:size(group, 2)+size(concatenated_array, 2)) = num2cell(concatenated_array);

rmData = cell2table(combined_cell,'VariableNames', variable_names_cell');
% Perform the repeated measures ANOVA
rm = fitrm(rmData, sprintf('s1-s%d ~ Group', day), 'WithinDesign', Sessions);
ranovaResults = ranova(rm);
fprintf('Repeated Measures ANOVA Results day %d CR percentage %s %s:\n', day, Gene, loc);
disp(ranovaResults);
end