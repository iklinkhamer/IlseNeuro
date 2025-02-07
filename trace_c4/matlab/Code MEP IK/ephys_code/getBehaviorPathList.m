%% IK 8-4-24
function path_list = getBehaviorPathList(mouse)

dir_content = dir(fullfile("/home/i.klinkhamer/Documents/Data/behaviorData/", mouse));
try
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    path_list = [];
end
if ~isempty(path_list)
    return
end

dir_content = dir(fullfile("/home/i.klinkhamer/Documents/Data/behaviorDataResults/", mouse));
try
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    path_list = [];
end
if ~isempty(path_list)
    return
end

trainingPath = "/home/i.klinkhamer/Documents/Data/behaviorDataTrainingOnlyMice/";
dir_content = dir(fullfile(trainingPath, mouse));
try
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    path_list = [];
end
if ~isempty(path_list)
    return

end
trainingPath = "/home/i.klinkhamer/Documents/Data/behaviorDataResultsTrainingOnlyMice/";
dir_content = dir(fullfile(trainingPath, mouse));
try
    path_list = IkUtils.getPathFromDir(dir_content);
catch
    path_list = [];
end
if ~isempty(path_list)
    return
else
    path_list = [];
end

end
