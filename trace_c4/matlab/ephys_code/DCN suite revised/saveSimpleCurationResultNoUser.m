%% I.K. 4-6-24
function ok = saveSimpleCurationResulNoUser(mname, masks)
    arguments
        mname(1,1) string
        masks(1,:) struct
    end
    
    simpleMasks = {masks.simple};
    saveSimpleMasks(mname, simpleMasks);

    modulationMasks = {masks.modulation};
    saveSimpleChannelModulationMasks(mname, modulationMasks);     
    
end