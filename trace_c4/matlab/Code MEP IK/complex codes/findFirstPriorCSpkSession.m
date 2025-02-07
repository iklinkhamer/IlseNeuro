function [counts, idcs] = findFirstPriorCSpkSession(mname)
    arguments
        mname (1,1) string = "Iowa";
    end
    
    sessionIdcs = splitSessionsByType(mname);
        
    if isempty(sessionIdcs.uniform)
        error ...
            ( "findFirstPriorCspkSession is not defined for %s because it has " ...
            + "no uniform sessions." ...
            , mname ...
            )
    end
    
    modulationMasks_ = loadComplexChannelModulationMasks(mname);
    modulationMasks = modulationMasks_(sessionIdcs.uniform);
    
    [counts, idcs] = cellfun(@countSession, modulationMasks, uniformoutput = false);
    
    ax = IkUtils.initPlots();
    plot(ax, 1:numel([counts{:}]), [counts{:}], '-*')
    title(ax, mname)
    xlabel(ax, 'session idx')
    ylabel(ax, 'prior-modulating CSpk count')
    
%     i = IkUtils.lazyFind(@checkSession, modulationMasks, onlyFirst = true);
    
    function [count, idcs] = countSession(channelMasks)
        
        if isempty(channelMasks)
            count = 0;
            idcs = [];
            return
        end
        mask = cellfun(@checkChannel, channelMasks);
        count = sum(mask);
        idcs = find(mask);

    end

%     function bool = checkSession(channelMasks)
%         
%         idx = IkUtils.lazyFind(@checkChannel, channelMasks, onlyFirst = true);
%         
%         bool = ~isempty(idx);
%         
%     end
    
    function bool = checkChannel(channelMask)
        
        if isstruct(channelMask) && isfield(channelMask, "prior")
            bool = channelMask.prior;
        else
            bool = false;
        end
        
    end
    
end