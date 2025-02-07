function filePath = getFilePath(mcode)
% P = IkUtils.getParams(); % IK change
% directory = P.dirHome;
% filePath = fullfile(directory, mname);
% addPaths(mfilename)
% structFilename = "StructEphysData.mat";
% stamp = nameDateStampFiles(mcode = mcode, s = 1, file_pattern = "StructEphysData", file_extension = ".mat");
[~, filepath_] = recursiveFileSearch("StructEphysData*.mat", mcode);
filePath = fileparts(filepath_(1));

% try
%     filename_try = fullfile(mcode, P.s(1), structFilename);
%     filePath_ = fileparts(which(filename_try));
%     filePath = erase(filePath_, P.s(1));
% catch
%     disp("Save location not found")
%     filePath = "";
%     return
% end

end