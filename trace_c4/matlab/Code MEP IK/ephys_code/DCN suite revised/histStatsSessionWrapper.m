%% I.K. 1-6-24
function stats = histStatsSessionWrapper ...
                            ( sessionData ...
                            , curation...
                            , kwargs ...
                            )
    arguments
        sessionData(1,1) struct
        curation(1,1) struct
        kwargs.event (1,1) string
        kwargs.modulating (1,1) logical
        kwargs.mname (1,1) string
    end
                
    if isempty(curation.simple) %isempty(curation.complex)
        stats = struct ...
            ( stats = [] ...
             , event = kwargs.event ...
             , modulating = kwargs.modulating ...
             , mname = kwargs.mname ...
            );
        return
    end
    
    p = IkUtils.getParams();

    function stats = computeStats(rasterField, roi, fullRange, mask)
        if not(any(mask))
            stats = struct.empty;
        else
            stats = arrayfun ...
            ( @(chanIdx) ...
                computePsthStats ...
                    ( sessionData.neuron(chanIdx).(rasterField) ...
                    , roi ...
                    , fullRange ...
                    ) ...
            , find(mask) ...
            );
        end
    end
    
    switch kwargs.event %lower(kwargs.event)
        case "cs_facilitation"
            rasterField = "RasterXY_cs_filtered";
            roi = p.sspkRanges.cs; %p.cspkRanges.cs;
            fullRange = p.psthRanges.cs_full;
            modulationMask = curation.modulation.cs_facilitation(:); %curation.modulation.cs(:);

        case "cs_suppression"
            rasterField = "RasterXY_cs_filtered";
            roi = p.sspkRanges.cs; %p.cspkRanges.cs;
            fullRange = p.psthRanges.cs_full;
            modulationMask = curation.modulation.cs_suppression(:); %curation.modulation.cs(:);
            
%         case "prior"
%             rasterField = "RasterXY_cs";
%             roi = p.cspkRanges.prior;
%             fullRange = p.psthRanges.cs_full;
%             modulationMask = curation.modulation.prior(:);
                        
        case "us_facilitation"
            rasterField = "RasterXY_us_filterd";
            roi = p.sspkRanges.us;
            fullRange = p.psthRanges.us_full;
            modulationMask = curation.modulation.us_facilitation(:);

        case "us_suppression"
            rasterField = "RasterXY_us_filterd";
            roi = p.sspkRanges.us;
            fullRange = p.psthRanges.us_full;
            modulationMask = curation.modulation.us_suppression(:);
            
        otherwise
            error("Unknown event: %s", kwargs.event)
    end
    
    if kwargs.modulating
        channelMask = curation.simple(:) & modulationMask; %curation.complex(:) & modulationMask;
    else
        channelMask = curation.simple(:) & not(modulationMask); %curation.complex(:) & not(modulationMask);
    end
    
    stats = struct ...
        ( stats = computeStats ...
                    ( rasterField ...
                    , roi ...
                    , fullRange ...
                    , channelMask ...
                    ) ...
        , event = kwargs.event ...
        , modulating = kwargs.modulating ...
        , mname = kwargs.mname ...
        );
end