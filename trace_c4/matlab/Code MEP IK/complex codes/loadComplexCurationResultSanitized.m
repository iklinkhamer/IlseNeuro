function [curationMasks, sessionsMask] = loadComplexCurationResultSanitized(mname, kwargs)
    arguments
        mname (1,1) string
        kwargs.onlyLatest (1,1) {islogical} = true;
    end    
   
    rawCurationResults = ...
        loadComplexCurationResultRaw(mname, onlyLatest = kwargs.onlyLatest);
    
    [curationMasks, sessionsMask] = arrayfun ...
        ( @sanitizeSession ...
        , rawCurationResults ...
        );
    
end

function [sessionMasks, isEmpty] = sanitizeSession(rawResult)
   
    comp = rawResult.complex;
    modul = rawResult.modulation;
    
    if isempty(comp)
        modulation = struct.empty;
    else
        modulation = struct ...
            ( cs = mkModulationMask(comp, modul, "cs") ...
            , prior =  mkModulationMask(comp, modul, "prior") ...
            , us = mkModulationMask(comp, modul, "us") ...
            );
    end
    
    sessionMasks = struct ...
        ( complex = comp ...
        , exemplary = rawResult.exemplary ...
        , modulation = modulation ...
        );
        
    isEmpty = isempty(comp);

end

function mask = mkModulationMask(channelMask, modulationCells, field)
    
    mask = false(size(channelMask));
    complexChannels = find(channelMask);
    
    for c = complexChannels(:)'
        if isempty(modulationCells{c})
            continue
        end
        mask(c) = modulationCells{c}.(field);
    end
    
end