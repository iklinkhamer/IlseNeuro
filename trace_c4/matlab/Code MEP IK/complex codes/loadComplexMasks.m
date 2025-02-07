function [masks, dates] = loadComplexMasks(mname, kwargs)
    arguments
        mname (1,1) string
        kwargs.onlyLatest (1,1) {islogical} = true;
    end     
   
    
    filePattern = sprintf ...
        ( "complexChannelMasks_%s" ...
        , mname ...
        );
    
    loadFn = @() IkUtils.io.loadTimestampedFiles ...
        ( filePattern ...
        , onlyLatest = kwargs.onlyLatest ...
        , folder = cspkResultsFolder() ...
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