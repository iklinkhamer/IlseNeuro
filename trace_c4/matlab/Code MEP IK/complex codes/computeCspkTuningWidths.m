function SDs = computeCspkTuningWidths(mname, kwargs)
    arguments
        mname (1,1) string
        kwargs.type (1,1) string = getMouseType(mname);
        kwargs.modulating (1,1) logical
        kwargs.event (1,1) string
    end
    
    import IkUtils.ifelse;
    
    switch lower(kwargs.type)
        case "uniform"
            sessionIdcs = splitSessionsByType(mname).uniform;
        case "delta"
            sessionIdcs = splitSessionsByType(mname).delta;
        case "deun"
            SDs = ...
                [ computeCspkTuningWidths( ...
                    mname ...
                    , type = "delta" ...
                    , modulating = kwargs.modulating ...
                    , event = kwargs.event ...
                    ) ...
                , computeCspkTuningWidths( ...
                    mname ...
                    , type = "uniform" ...
                    , modulating = kwargs.modulating ...
                    , event = kwargs.event ...
                    ) ...
                ];
            return
        otherwise
            error("Unknown session type: %s", kwargs.type)
    end
    
    allStats = histStatsMouseWrapper ...
        ( mname ...
        , event = kwargs.event ...
        , modulating = kwargs.modulating ...
        );
    statsData = allStats(sessionIdcs);

    eventData = IkUtils.flatten({statsData.stats}, uniformOutput = true);
    
    if isempty(eventData)
        widths = [];
    else    
        [~, filterMask] = rejectSDOutliers([eventData.maxAmpTime] * 1000);
        widths = [eventData(filterMask).sd] * 1000;
    end
    
%     fprintf ...
%         ( "Inspecting %s %s, %s %s, found %d channels of which %d outliers.\n" ...
%         , mname ...
%         , kwargs.type ...
%         , kwargs.event ...
%         , ifelse(kwargs.modulating, "modulating", "non modulating") ...
%         , numel(eventData) ...
%         , numel(eventData) - numel(widths) ...
%         );
    
    SDs = struct ...
        ( sds = widths ...
        , sds_mean = mean(widths(~isnan(widths))) ...
        , type = kwargs.type ...
        , event = kwargs.event ...
        , modulating = kwargs.modulating ...
        , mname = mname ...
        );
    
end