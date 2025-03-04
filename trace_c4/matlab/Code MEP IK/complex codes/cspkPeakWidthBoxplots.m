function [res, axs, fig] = cspkPeakWidthBoxplots(kwargs)
    arguments
        kwargs.event (1,1) string = "cs";
        kwargs.modulating (1,1) logical = true;
        kwargs.mouseNames (1,:) string = [defaultMice().name];
    end
    
    import IkUtils.flatmap;
    import IkUtils.ifelse;
    
%     mice = defaultMice();
    mouseNames = kwargs.mouseNames; % [mice.name];
    peakWidths = flatmap ...
        ( @(mname) computeCspkTuningWidths ...
                    ( mname ...
                    , modulating = kwargs.modulating ...
                    , event = kwargs.event ...
                    ) ...
        , mouseNames ...
        );
    
    deltaMask = strcmpi([peakWidths.type], "delta");
    uniformMask = not(deltaMask);
    
    [axs, fig] = IkUtils.initPlots([1 5], layout = {{1 1 1 3} {1 4 1 2}});
    
    ylabel(axs(1), "Mean of peak widths across neurons")
        
    deltaData = [peakWidths(deltaMask).sds_mean];
    deltaSubjects = [peakWidths(deltaMask).mname];
    uniformData = [peakWidths(uniformMask).sds_mean];
    uniformSubjects = [peakWidths(uniformMask).mname];
    
    [h, p, ci, stats] = ttest2(deltaData, uniformData, tail = "left")
    
    tTestResult = struct ...
        ( h = h ...
        , p = p ...2
        , ci = ci ...
        , stats = stats ...
        );
    
    mkDeltaVsUniformBoxPlot ...
        ( axs ...
        , deltaData ...
        , uniformData ...
        , sprintf ...
            ( "%s - %s" ...
            , kwargs.event ...
            , ifelse(kwargs.modulating, "modulating", "non modulating") ...
            ) ...
        , deltaSubjects = deltaSubjects ...
        , uniformSubjects = uniformSubjects ...
        , tTestResult = tTestResult ...
        )
    
    res = struct ...
        ( delta = deltaData ...
        , uniform = uniformData ...
        , event  = kwargs.event ...
        , modulating = kwargs.modulating ...
        , ttest = tTestResult ...
        );
    
end