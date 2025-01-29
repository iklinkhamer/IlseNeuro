%% IK 8-4-24
function plotErrorBarFigure(Data, ID, fig, subcoords, kwargs)
arguments
    Data = [];
    ID = [];
    fig = [];
    subcoords = [1 1 1];
    kwargs.showBackGroundLines = true;
    kwargs.day = min(IkUtils.getParams().lastDay, sum(any(~isnan(Data.onsets), 2)));
end
% if nargin < 4
%     showBackGroundLines = 0;
% end
showBackGroundLines = kwargs.showBackGroundLines;
day = kwargs.day;
if length(Data) < 2
%     lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(Data.onsets), 2)));
%     if nargin < 5
%         day = lastDay;
%     end
    axes(fig); hold on
    % set(gcf, 'DefaultFigureRenderer', 'painters');
    subplot(subcoords(1), subcoords(2), subcoords(3))
    hold on
    errorbar(Data.crampmean,Data.crampSEM, 'b','LineWidth',1.2);
    xticks(1:day);
    xlim([1, day]);
    ylim([0 100]);
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR amplitude')
    title('CR amplitude during training')
    box('off');
    hold off

else
    WT = Data(1);
    MUT = Data(2);
%     lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(WT.onsets), 2)));

%     if nargin < 5
%         day = lastDay;
%     end

    axes(fig); hold on
    % set(gcf, 'DefaultFigureRenderer', 'painters');
    subplot(subcoords(1), subcoords(2), subcoords(3))
    hold on
    if showBackGroundLines
        plot(nanmean(WT.cramp, 3), 'color', [0.68, 0.85, 0.9, 0.5]);
        plot(nanmean(MUT.cramp, 3), 'color', [1, 0.8,0.8, 0.5]);
    end
    e1 = errorbar(WT.crampmean,WT.crampSEM, 'b','LineWidth',1.2);
    e2 = errorbar(MUT.crampmean,MUT.crampSEM, 'r','LineWidth',1.2);
    xticks(1:day);
    xlim([1, day]);
    ylim([0 100]);
    set(gca,'TickDir','out');
    xlabel('Session');
    ylabel('CR amplitude')

    if length(ID) > 1
        Gene = ID(1);
        loc = ID(2);
        title(sprintf('CR amplitude %s %s', Gene, loc))
    else
        title(sprintf('CR amplitude %s', ID))
    end

    box('off');
    legend([e1, e2],{sprintf('WT %s', loc), sprintf('MUT %s', loc)})
    hold off

    repmAnova(WT,MUT,day, Gene, loc)
end
end

function repmAnova(WT,MUT, day, Gene, loc)
P = IkUtils.getParams();
CSmin = P.tracesRanges.cs.min;
CSmax = P.tracesRanges.cs.max;
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
Sessions = table((1:day)', 'VariableNames', {'Sessions'});

WT_cramp = WT.cramp(1:day,:,:);
MUT_cramp = MUT.cramp(1:day,:,:);

WT_cramp_mean_mice = nanmean(WT_cramp,3);
MUT_cramp_mean_mice = nanmean(MUT_cramp,3);

[~,colWT] = find(isnan(WT_cramp_mean_mice));
WT_cramp(:,colWT,:) = [];
[~,colMUT] = find(isnan(MUT_cramp_mean_mice));
MUT_cramp(:,colMUT,:) = [];

variable_names = [{'Group'}, strcat('s', string(1:day))];
variable_names_cell = cellfun(@(x) {x}, variable_names);


for s = 1:size(WT_cramp,1)
    for m = 1:size(WT_cramp,2)
        traces_session = WT.traces(s,(m-1)*P.nTrials + 1 : m*P.nTrials, :);
        baseline_min_CS = nanmin(nanmean(traces_session(:,WT.type(m,:) ==0,CSmin:CSmax),3));
        for t = 1:size(WT_cramp, 3)
            if isnan(WT_cramp(s,m,t))
                WT_cramp(s,m,t) = baseline_min_CS;
            end
        end
    end
end

for s = 1:size(MUT_cramp,1)
    for m = 1:size(MUT_cramp,2)
        traces_session = MUT.traces(s,(m-1)*P.nTrials + 1 : m*P.nTrials, :);
        baseline_min_CS = nanmin(nanmean(traces_session(:,MUT.type(m,:) ==0,CSmin:CSmax),3));
        for t = 1:size(MUT_cramp, 3)
            if isnan(MUT_cramp(s,m,t))
                MUT_cramp(s,m,t) = baseline_min_CS;
            end
        end
    end
end

[sessionsWT, miceWT, trialsWT] = size(WT_cramp);
reshaped_data_WT_new = reshape(WT_cramp, [sessionsWT, miceWT*trialsWT]);
groupWT = repmat({'WT'}, size(reshaped_data_WT_new, 2),1);

[sessionsMUT, miceMUT, trialsMUT] = size(MUT_cramp);
reshaped_data_MUT = reshape(MUT_cramp, [sessionsMUT, miceMUT*trialsMUT]);
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
fprintf('Repeated Measures ANOVA Results day %d CR amplitude %s %s:\n', day, Gene, loc);
disp(ranovaResults);
end

