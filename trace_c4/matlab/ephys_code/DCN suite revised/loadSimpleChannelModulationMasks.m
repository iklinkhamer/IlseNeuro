%% I.K. 1-6-24
function [masks, dates] = loadSimpleChannelModulationMasks(mcode, kwargs)
    arguments
        mcode (1,1) string
        kwargs.onlyLatest (1,1) {islogical} = true;
    end     
   
    
    prefixPattern = sprintf ...
        ( "simpleChannelModulationMasks_%s" ...
        , mcode ...
        );

    folder_ = sspkResultsFolder(prefixPattern, IkUtils.getParams().pathSpikeSortingHome);
    folder = erase(folder_, mcode);

    loadFn = @() IkUtils.io.loadTimestampedFiles ...
        ( prefixPattern ...
        , mcode... 
        , onlyLatest = kwargs.onlyLatest ...  
        , folder = folder ...
        , uniformData = false ... % channel masks are of varying length
        );
    
    if nargout < 2
        % Don't compute date strings if they weren't requested by the
        % caller
        masks_ = loadFn();
    else
        [masks_, dates] = loadFn();
    end
    
    if isempty(masks_)
        masks = masks_;
        return
    end
    
    masks = masks_{1}.masks;
    
end