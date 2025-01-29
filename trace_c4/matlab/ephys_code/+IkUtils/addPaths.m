%% IK 8-4-24
function addPaths(currentFileName)
if nargin < 1
    currentFileName = mfilename;
end
path_current = fileparts(which(currentFileName));
addpath(genpath(path_current));

% Add Data folders
p = IkUtils.getParams;
addpath(genpath(p.dirHomeData));
addpath(genpath(p.dirDATAData));
addpath(genpath(p.dirHomeCode));
end