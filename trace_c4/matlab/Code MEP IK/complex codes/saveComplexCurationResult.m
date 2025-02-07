% TODO: check for existence of the requested files on disk after saving and
% report this as a logical return value
function ok = saveComplexCurationResult(mname, masks)
    arguments
        mname(1,1) string
        masks(1,:) struct
    end
    
    complexMasks = {masks.complex};
    
    saveComplexMasks(mname, complexMasks);
    
    exemplaryMasks = {masks.exemplary};
    saveExemplaryComplexChannelMasks(mname, exemplaryMasks);
    
    modulationMasks = {masks.modulation};
    saveComplexChannelModulationMasks(mname, modulationMasks);
    
    
    
end