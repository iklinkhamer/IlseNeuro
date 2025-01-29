%% I.K. 1-6-24
function [masks, dates] = loadTimestampedFiles(prefixPattern, mcode, kwargs)
arguments
    prefixPattern = "";
    mcode = "";
    kwargs.onlyLatest = true;
    kwargs.folder = IkUtils.getParams().dirHomeData;
    kwargs.uniformData = false;
end
dir_contents = dir(kwargs.folder);
if ~isempty(dir_contents)
    path_list = IkUtils.getPathFromDir(dir_contents);
else
    masks = [];
    dates = [];
    return
end

for f = path_list
    if contains(f, mcode)
        filename = fullfile(f, sprintf("%s*",prefixPattern));
    else
        filename = fullfile(f, mcode, sprintf("%s*",prefixPattern));
    end
    dir_files = dir(filename);
    if ~isempty(dir_files)
        dir_contents_file = dir_files;
        break
    end
end

if ~isempty(dir_contents_file)
    if kwargs.onlyLatest
        try
            masks = {load(fullfile(kwargs.folder, mcode, dir_contents_file(end).name))}; 
            dates = erase(dir_contents_file(end).name, prefixPattern);
        catch
            masks = [];
            dates = [];
        end
    else
        tic
        for i = 1:length(dir_contents_file)
            masks(i) = load(fullfile(kwargs.folder, mcode, dir_contents_file(i).name));
            dates = erase(dir_contents_file(i).name, prefixPattern);
            if toc > 3
                continue
            end
        end
    end
else
    masks = [];
    dates = [];
end
end
