%% IK 25-4-24       Recursively look for files containing a certain pattern in their filename.
function [files, folders] = recursiveFileSearch(filename_pattern, mcode, s, parent_directory, max_depth)
arguments
    filename_pattern = ''; %"Config_h32_oe_*.csv";
    mcode = ''; %"20-MI19442-05"
    s = []; %1 % You can enter either a single session or an array of sessions. 
    parent_directory = [IkUtils.getParams().dirHomeData, IkUtils.getParams().dirDATAData];
    max_depth = 3;
end
P = IkUtils.getParams();
matching_files = [];
if isempty(filename_pattern)
    disp("Didn't enter filename pattern")
    files = [];
    folders = [];
    return
end
files = [];
folders = [];
for d = 1 : length(parent_directory)
    pDir = parent_directory(d);
    directories = dir(pDir);
    files_dir = directories(~[directories.isdir]);
    directories = directories([directories.isdir] & ~strcmp({directories.name}, ".") & ~strcmp({directories.name}, ".."));

    matching_files = searchFilesRecursively(matching_files, directories, filename_pattern, max_depth);  % Recursively search for files matching the pattern in all subdirectories

    partsFilenamePattern = split(filename_pattern, "*");
    for f = 1 : numel(files_dir)
        if all(arrayfun(@(p) contains(files_dir(f).name, p), partsFilenamePattern))
            matching_files = [matching_files, files_dir(f)];
        end
    end

    if ~isempty(matching_files)
        if ~isempty(mcode)
            mcodeMask = contains(string({matching_files.name}), mcode) | contains(string({matching_files.folder}), mcode);
            if any(mcodeMask)
                matching_files = matching_files(mcodeMask);
            else
                continue
            end
        end

        if ~isempty(s)
            if all(arrayfun(@(x) isfloat(x) | isinteger(x), s)) && all(s <= P.n_sessions)
                sFolderMask = arrayfun(@(session) contains(string({matching_files.name}), P.s(session)) | contains(string({matching_files.folder}), P.s(session)), s, 'UniformOutput', false);
                sFolderMask = any(cat(3, sFolderMask{:}), 3);
                if ~isempty(mcode) && sum(sFolderMask) == 0
                    mice = defaultMice();
                    ephysDates = mice([defaultMice().code] == mcode).ephysdates;
                    if ~isempty(ephysDates)
                        ephysDates = ephysDates(s);
                        sEphysMask = arrayfun(@(date) contains(string({matching_files.name}), date) | contains(string({matching_files.folder}), date), ephysDates, 'UniformOutput', false);
                        sEphysMask = any(cat(3, sEphysMask{:}), 3);

                        parts = split(ephysDates, "-"); 
                        if numel(s) > 1
                            ephysDatesJoined = join(parts, "", 3);
                        else
                            ephysDatesJoined = join(parts, "");
                        end
                        ephysDatesJoined(ephysDatesJoined.strlength > 6) = extractAfter(ephysDatesJoined(ephysDatesJoined.strlength > 6), 2);
                        sJoinedMask = arrayfun(@(date) contains(string({matching_files.name}), date) | contains(string({matching_files.folder}), date), ephysDatesJoined, 'UniformOutput', false);
                        sJoinedMask = any(cat(3, sJoinedMask{:}), 3);
                        fullDateMask = sEphysMask | sJoinedMask;     
                    else
                        fullDateMask = sFolderMask;
                    end
                else
                    fullDateMask = sFolderMask;
                end
            else
                s = string(s); % Convert into string in case it wasn't already a string.
                s(s.strlength > 7) = extractAfter(s(s.strlength > 7), 2);
                sStringMask = arrayfun(@(date) contains(string({matching_files.name}), date) | contains(string({matching_files.folder}), date), s, 'UniformOutput', false);
                sStringMask = any(cat(3, sStringMask{:}), 3);

                fullDateMask = sStringMask;
            end
            if any(fullDateMask)
                matching_files = matching_files(fullDateMask);
            else
                continue
            end
        end

        files = string({matching_files.name});
        folders = string({matching_files.folder});
    end
end

end

function [matching_files, directories_new] = searchFiles(matching_files, directory, filename_pattern)
folder_path = fullfile(directory.folder, directory.name);
matching_files_in_folder = dir(fullfile(folder_path, filename_pattern));
matching_files = [matching_files; matching_files_in_folder];
directories_new = dir(folder_path);
directories_new = directories_new([directories_new.isdir] & ~strcmp({directories_new.name}, ".") & ~strcmp({directories_new.name}, ".."));
end

function matching_files = searchFilesRecursively(matching_files, directories, filename_pattern, depth)
if depth == 0
    return;
end

for i = 1:numel(directories)
    [matching_files, next_directories] = searchFiles(matching_files, directories(i), filename_pattern);
    matching_files = searchFilesRecursively(matching_files, next_directories, filename_pattern, depth - 1);
end
end