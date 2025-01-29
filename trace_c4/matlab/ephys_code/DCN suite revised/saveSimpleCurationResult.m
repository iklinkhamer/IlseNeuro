%% I.K. 1-6-24
function ok = saveSimpleCurationResult(mname, masks)
    arguments
        mname(1,1) string
        masks(1,:) struct
    end
    
    simpleMasks = {masks.simple};
    saveSimpleMasks(mname, simpleMasks);
    
    exemplaryMasks = {masks.exemplary};
    saveExemplarySimpleChannelMasks(mname, exemplaryMasks);
    
    modulationMasks = {masks.modulation};
    saveSimpleChannelModulationMasks(mname, modulationMasks);   
    
    
end