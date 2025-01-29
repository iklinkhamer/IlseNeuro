%% IK 5-4-24
function smoothedTraces = smootheTraces(traces)
% Parameters
windowSize = 7; % Window size for Savitzky-Golay filter
degree = 2; % Degree of polynomial for Savitzky-Golay filter

smoothedTraces = sgolayfilt(traces', degree, windowSize);

alpha = 0.05; % Significance level for Grubb's test
maxIterations = 5; % Maximum number of iterations for outlier removal
% Iterative Grubb's outlier detection
baselineSDs = std(smoothedTraces(1:40, :));
baselineMeans = mean(smoothedTraces(1:40, :));
isOutlier = false(size(baselineSDs));
for iter = 1:maxIterations
    meanSD = mean(baselineSDs);
    stdSD = std(baselineSDs);
    zScoresSD = abs((baselineSDs - meanSD) / stdSD);
    [~, idxSD] = max(zScoresSD);

    meanmean = mean(baselineMeans);
    stdMean = std(baselineMeans);
    zScoresMean = abs((baselineMeans - meanmean) / stdMean);
    [~, idxMean] = max(zScoresMean);

    zScoreThreshold = tinv(1 - alpha / (2 * numel(baselineSDs)), numel(baselineSDs) - 1);

    if zScoresSD(idxSD) > zScoreThreshold
        isOutlier(idxSD) = true;
        baselineSDs(idxSD) = [];
        baselineMeans(idxSD) = [];
        if idxSD < idxMean
            idxMean = idxMean - 1;
        end
    end
    if zScoresMean(idxMean) > zScoreThreshold
        isOutlier(idxMean) = true;
        baselineSDs(idxMean) = [];
        baselineMeans(idxMean) = [];
    else
        break;
    end
end
% Remove trials with unstable baselines
if any(isOutlier)
    smoothedTraces(:, isOutlier) = nan;
end
smoothedTraces = smoothedTraces';
end