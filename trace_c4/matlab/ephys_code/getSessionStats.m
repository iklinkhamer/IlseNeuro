function sessionStats = getSessionStats(curationResults)
histogramStats = [curationResults.histogramStats];

distTypesSessions = [curationResults.distType];
histColor = getHistColor(distTypesSessions(1));

sessionStats = struct;

for s = 1 : length(histogramStats)
    modulatingIdcs = find(curationResults(s).modulatingMask);
    try
        nonModulatingIdcs = find...
            ( curationResults(s).modulatingMask(:) == 0 ...
            & curationResults(s).complexMasks(:) == 1 ...
            );
    catch
    end
    try
        nonModulatingIdcs = find...
            ( curationResults(s).modulatingMask(:) == 0 ...
            & curationResults(s).simpleMasks(:) == 1 ...
            );
    catch
    end
%     nonModulatingIdcsMask = ones(length(histogramStats(s).max_amp),1);
%     nonModulatingIdcsMask(modulatingIdcs) = 0;
%     nonModulatingIdcs = find(nonModulatingIdcsMask == 1);
    
    sessionStats(s).peakTimesMod = histogramStats(s).max_amp_delay_time(modulatingIdcs);
    sessionStats(s).peakTimesNonMod = histogramStats(s).max_amp_delay_time(nonModulatingIdcs);
    sessionStats(s).sdMod = histogramStats(s).st_dev(modulatingIdcs);
    sessionStats(s).sdNonMod = histogramStats(s).st_dev(nonModulatingIdcs);
end
end