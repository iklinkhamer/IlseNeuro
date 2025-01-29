%% I.K. 5-6-24
% Called by function: plotSummaryFiguresMUTvsWT_IK.m
function [boxplot_data, boxplot_data_cell] = getBoxplotData(Data, field)

if length(Data) > 1
    boxplot_data = [];
    boxplot_data_cell = [];
    stimuli = ["cs", "us"];
    for s = 1:length(stimuli)
        stimulus = stimuli(s);
        allStatsWT = Data(1).(stimulus);
        allStatsMUT = Data(2).(stimulus);

        facFieldName = stimulus + "_facilitation";
        supFieldName = stimulus + "_suppression";

        maskFacilitationMUT = [allStatsMUT.(facFieldName)] == 1;
        maskFacilitationWT = [allStatsWT.(facFieldName)] == 1;
        maskSuppressionMUT = [allStatsMUT.(supFieldName)] == 1;
        maskSuppressionWT = [allStatsWT.(supFieldName)] == 1;
        maskNonModMUT = [allStatsMUT.(facFieldName)] == 0 & [allStatsMUT.(supFieldName)] == 0;
        maskNonModWT = [allStatsWT.(facFieldName)] == 0 & [allStatsWT.(supFieldName)] == 0;

        boxplot_data_stim = [ ...
            [allStatsMUT(maskNonModMUT).(field)] ...
            [allStatsWT(maskNonModWT).(field)] ...
            [allStatsMUT(maskFacilitationMUT).(field)] ...
            [allStatsWT(maskFacilitationWT).(field)] ...
            [allStatsMUT(maskSuppressionMUT).(field)] ...
            [allStatsWT(maskSuppressionWT).(field)] ...
            ];
        boxplot_data_cell_stim = { ...
            [allStatsMUT(maskNonModMUT).(field)]; ...
            [allStatsWT(maskNonModWT).(field)]; ...
            [allStatsMUT(maskFacilitationMUT).(field)]; ...
            [allStatsWT(maskFacilitationWT).(field)]; ...
            [allStatsMUT(maskSuppressionMUT).(field)]; ...
            [allStatsWT(maskSuppressionWT).(field)] ...
            };

        boxplot_data = [boxplot_data; {boxplot_data_stim}];
        boxplot_data_cell = [boxplot_data_cell, boxplot_data_cell_stim];
    end

else
    boxplot_data = [];
    boxplot_data_cell = [];
    stimuli = ["cs", "us"];
    for s = 1:length(stimuli)
        stimulus = stimuli(s);
        allStats = Data.(stimulus);

        facFieldName = stimulus + "_facilitation";
        supFieldName = stimulus + "_suppression";

        maskFacilitationWT = [allStats.(facFieldName)] == 1;
        maskSuppressionWT = [allStats.(supFieldName)] == 1;
        maskNonModWT = [allStats.(facFieldName)] == 0 & [allStats.(supFieldName)] == 0;

        boxplot_data_stim = [ ...
            [allStats(maskNonModWT).(field)] ...
            [allStats(maskFacilitationWT).(field)] ...
            [allStats(maskSuppressionWT).(field)] ...
            ];
        boxplot_data_cell_stim = { ...
            [allStats(maskNonModWT).(field)]; ...
            [allStats(maskFacilitationWT).(field)]; ...
            [allStats(maskSuppressionWT).(field)] ...
            };
        boxplot_data = [boxplot_data; {boxplot_data_stim}];
        boxplot_data_cell = [boxplot_data_cell, boxplot_data_cell_stim];
    end

end

end