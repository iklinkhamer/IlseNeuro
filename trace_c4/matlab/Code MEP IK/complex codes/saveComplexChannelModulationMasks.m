function saveModulationMasks(mname, masks)
    
    filePath = fullfile ...
        ( cspkResultsFolder() ...
        , sprintf ...
            ( "complexChannelModulationMasks_%s_%d.mat" ...
            , mname ...
            , IkUtils.now ...
            ) ...
        );
        
    save(filePath, "masks")
    
end