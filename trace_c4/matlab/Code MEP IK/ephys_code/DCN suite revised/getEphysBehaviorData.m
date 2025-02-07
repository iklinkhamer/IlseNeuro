% TODO: This function contains unnecessary specialization. Unify and centralize
% this with the other data loading code.
function behaviorData = getEphysBehaviorData(fileList)
    
    behaviorData = arrayfun ...
        ( @(file) load(fullfile(file.folder,file.name)) ...
        , fileList ...
        );
    
end