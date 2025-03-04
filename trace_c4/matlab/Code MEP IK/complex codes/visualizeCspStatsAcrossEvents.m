function visualizeCspStatsAcrossEvents(mouseNames)
    arguments
        mouseNames (1,:) string = [paperFiguresMice().name];
    end
    
    uiP = IkUtils.visualization.Params();
    p = getParams();
    events = ...
        [ "cs" ... % IK removed prior
        , "us" ...
        ];
    eventDelayOffsetCorrections = ...
        [ 0 ...
        , p.delayEventsCs.Uniform(2) * 1000 ...
        , 0 ...
        ];
    
    [rangeMins, rangeMaxes] = structfun ...
        ( @(event) deal(event.min * 1000, event.max * 1000) ...
        , p.cspkRanges ...
        );
    rangeMin = min(rangeMins);
    rangeMax = max(rangeMaxes - rangeMins);
    
    
        
    delayDataWrapper = @(event, correction) unpack ...
        ( arrayfun ...
            ( @(mname) computeCspkTuningDelays ...
                ( mname ...
                , modulating = true ...
                , event = event ...
                , type = "uniform" ...
                ) ...
            , mouseNames ...
            ) ...
        , "delays" ...
        ) - correction;
    ...- p.cspkRanges.(event).min * 1000;
    
    delayDataCells = arrayfun ...
        ( delayDataWrapper ...
        , events ...
        , eventDelayOffsetCorrections ...
        , uniformoutput = false ...
        );
    
    
    widthDataWrapper = @(event) unpack ...
        ( arrayfun ...
            ( @(mname) computeCspkTuningWidths ...
                ( mname ...
                , modulating = true ...
                , event = event ...
                , type = "uniform" ...
                ) ...
            , mouseNames ...
            ) ...
        , "sds" ...
        );
    
    widthDataCells = arrayfun ...
        ( widthDataWrapper ...
        , events ...
        , uniformoutput = false ...
        );
    
    nChannelsTotal = sum( cellfun( @numel, widthDataCells ) );
    
    [delayAx, delayFig] = processDataSet(delayDataCells, "mode");
    [widthAx, widthFig] = processDataSet(widthDataCells, "SD");

    
    
    function [ax, fig] = processDataSet(dataCells, analysis)
    

        [ax, fig] = IkUtils.initPlots ...
            ( [10 1] ...
            , layout = {{2 1 9 1}} ...
            );
        
        rangeLengths = (rangeMaxes - rangeMins)';
    
        arrayfun ...
            ( @(idx, rngLen) ...
                plotSingleEvent ...
                    ( ax ...
                    , dataCells{idx} ...
                    , idx ...
                    , rngLen ...
                    ) ...
            , 1:3 ...
            , rangeLengths ...
            )
        
        title ...
            ( ax ...
            , sprintf("CSpk %s", analysis) ...
            )
        subtitle ...
            ( ax ...
            , sprintf("N=%d", nChannelsTotal) ...
            )
        ylabel(ax, sprintf("CSpk histogram %s [ms]", analysis))
        xlim(ax, [.5 3.5])
        ax.XTick = 1:numel(events);
        ax.XTickLabel = cellstr(events);
        
        plotEventPairTtest(dataCells{1:2}, [1,2])
        plotEventPairTtest(dataCells{2:3}, [2,3])

        function plotEventPairTtest(leftData, rightData, idcs)
        
            [h, pVal] = ttest2(leftData, rightData)
            
keyboard
            
            xLength = idcs(2) - idcs(1);
            xs = ...
                [ idcs(1) + xLength/20 ...
                , idcs(2) - xLength/20 ...
                ];
            y = max(rangeLengths(idcs)) + 15;
            
            l = line ...
                ( ax ...
                , xs, [y y] ...
                , color = uiP.foregroundColor ...
                , clipping = 'off' ...
                , handlevisibility = 'off' ...
                );
            
            if pVal > 5e-2
                str = "n.s.";
            elseif pVal <= 1e-2
                str = "**";
            else
                str = "*";
            end
            
            t = text ...
                ( ax ...
                , mean(xs), y+6 ...
                , str ...
                , color = uiP.foregroundColor ...
                , horizontalalignment = 'center' ...
                );
                        
        end
        
        function plotSingleEvent(ax, data, xPos, yMax)
keyboard
            avg = mean(data);
            stdErr = std(data); % / sqrt(numel(delays));
        
            bar ...
                ( ax ...
                , xPos ...
                , avg ...
                , facecolor = [.85 .85 .85] ...
                , edgecolor = [.7 .7 .7] ...
                )
            line ...
                ( ax ...
                , repmat(xPos, size(data)) + randn(size(data))/20 ...
                , data ...
                , marker = '+' ...
                , linestyle = 'none' ...
                , color = [.7 .7 .7] ...
                )
            errorbar ...
                ( ax ...
                , xPos ...
                , avg ...
                , stdErr ...
                , color = [.35 .35 .35] ...
                )
            line ...
                ( ax ...
                , xPos + [-0.3 0.3] ...
                , [yMax yMax] ...
                , linestyle = ':' ...
                , color = [.35 .35 .35] ...
                )
            
            text ...
                ( ax ...
                , xPos, yMax + 5 ...
                , sprintf("N=%d", numel(data)) ...
                , HorizontalAlignment = 'center' ...
                , color = [.3 .3 .3] ...
                );
        
        end
    
    end
    
end
    
function [values, labels] = unpack(data, dataField)
    
    import IkUtils.flatmap;

    labels = flatmap ...
        ( @(dat) repmat(dat.mname, size(dat.(dataField))) ...
        , data ...
        );
    values = [data.(dataField)];

end