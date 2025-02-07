%% I.K. 1-6-24
function [curationMasks, sessionsMask] = loadSimpleCurationResultSanitized(mname, kwargs)
    arguments
        mname (1,1) string
        kwargs.onlyLatest (1,1) {islogical} = true;
    end
    
    rawCurationResults = ...
        loadSimpleCurationResultRaw(mname, onlyLatest = kwargs.onlyLatest);
    
    [curationMasks, sessionsMask] = arrayfun ...
        ( @sanitizeSession ...
        , rawCurationResults ...
        );
    
end

function [sessionMasks, isEmpty] = sanitizeSession(rawResult)
   
    simp = rawResult.simple;
    modul = rawResult.modulation;
    
    sessionMasks = struct ...
        ( simple = simp ...
        , exemplary = rawResult.exemplary ...
        , modulation = struct ...
            ( cs_facilitation = mkModulationMask(simp, modul, "cs_facilitation") ... % IK removed prior
            , cs_suppression = mkModulationMask(simp, modul, "cs_suppression") ...      
            , us_facilitation = mkModulationMask(simp, modul, "us_facilitation") ...
            , us_suppression = mkModulationMask(simp, modul, "us_suppression") ...
            ) ...
        );
        
    isEmpty = isempty(simp);

%         , rawResult.simple ...
%         , rawResult.exemplary ...
%         , rawResult.modulation ...
    
end

function mask = mkModulationMask(channelMask, modulationCells, field)
    
    mask = false(size(channelMask));
    simpleChannels = find(channelMask);
    
    for c = simpleChannels(:)'
        if isempty(modulationCells{c})
            continue
        end
        mask(c) = modulationCells{c}.(field);
    end
    
end