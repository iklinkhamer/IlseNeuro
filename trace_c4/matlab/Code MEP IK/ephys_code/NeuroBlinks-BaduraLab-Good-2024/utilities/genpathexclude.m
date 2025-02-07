function generated_path = genpathexclude(basefolder, exclude_pattern)
% Like built-in genpath but allows you to exclude certain folders from generated path, specified as pattern

generated_path = genpath(basefolder);



generated_path = 