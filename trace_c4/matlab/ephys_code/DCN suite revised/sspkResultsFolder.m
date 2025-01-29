%% IK 8-4-24
function sspkFolder = sspkResultsFolder(prefixPattern, dirPath)
if nargin < 2
    dirPath = IkUtils.getParams().pathSpikeSortingHome;
end
directories = strsplit(genpath(dirPath), pathsep);
for i = 1:length(directories)
    current_directory = directories{i};
    dir_files = dir(fullfile(current_directory, sprintf("%s*", prefixPattern)));    
    if ~isempty(dir_files)
        sspkFolder = fileparts(dir_files(1).folder);
        return
    else
        continue
    end
end
if ~exist("sspkFolder", "var")
    sspkFolder = '';
end
        
end