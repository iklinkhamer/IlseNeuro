%% I.K. 1-6-24
function [masks, dates] = loadSimpleMasks(mcode, kwargs)
    arguments
        mcode (1,1) string
        kwargs.onlyLatest (1,1) {islogical} = true;
    end     
    
    prefixPattern = sprintf ...
        ( "simpleChannelMasks_%s" ...
        , mcode ...
        );
    
    loadFn = @() IkUtils.io.loadTimestampedFiles ...
        ( prefixPattern ...
        , mcode... % IK change
        , onlyLatest = kwargs.onlyLatest ...
        , folder = erase(sspkResultsFolder(prefixPattern), mcode) ...
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