%% IK 8-4-24
function path_list = getPathFromDir(dir_contents)

f = 1;
for i = 1:length(dir_contents)
    
    if dir_contents(i).name == "." || string(dir_contents(i).name) == ".."
        continue;
    end
    foldername = string(dir_contents(i).name);
    path = string(dir_contents(i).folder);

    path_list(f) = fullfile(path, foldername);
    f = f+1;
end

end
