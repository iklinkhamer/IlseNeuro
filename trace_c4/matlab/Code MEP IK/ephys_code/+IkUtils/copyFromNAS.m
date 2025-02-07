%% IK 5-4-24
function copyFromNAS()

app = uimain(scriptname=mfilename);
mouseCodes = [defaultMice().code];
mouseNames = [defaultMice().name];
mouseNames = [mouseNames, "N/A"];
mouseCodes = [mouseCodes, "N/A"];
app.MouseDropDown.Items = mouseNames;
app.StartmouseDropDown.Items = mouseNames;
app.CodeDropDown.Items = mouseCodes;
app.DropDown.Items = mouseCodes;

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    appData = evalin('base', 'appData');

    if appData.enumeration == "Single mouse"
        mousecodes = appData.mouseCode;
    elseif appData.enumeration == "All mice of a single group"
        Type = string(appData.gene) + string(appData.loc) + string(appData.type);
        mouseTypes = [defaultMice().type];
        ephys_mask = arrayfun(@(x) ~isempty(x.ephysdates), defaultMice());
        mousecodes = mouseCodes(mouseTypes == Type & ephys_mask);
        if contains(mousecodes, appData.startMouseCode)
            mousecodes = mousecodes(mousecodes == appData.startMouseCode:end);
        end
    elseif appData.enumeration == "All mice"
        ephys_mask = arrayfun(@(x) ~isempty(x.ephysdates), defaultMice());
        mousecodes = mouseCodes(ephys_mask);
        if contains(mousecodes, appData.startMouseCode)
            mousecodes = mousecodes(mousecodes == appData.startMouseCode:end);
        end
    end
    startFolder = appData.startFolder;

    mousecodes(mousecodes == "21-MI10159-01"|mousecodes == "21-MI10159-02"|mousecodes == "21-MI10159-06"|mousecodes == "Shank2KOMUT"|mousecodes == "Shank2KOWT") = [];

    sessionDates = struct;
    for m = 1 : length(mousecodes)
        mcode = mousecodes(m);
        sessionDates(m).mouse = mcode;
        fprintf("Starting copying for mouse: %s\n", mcode)
        mcodeParts = split(mcode, ".");
        if length(mcodeParts) > 1
            part1_ = char(mcodeParts(1));
            part1 = string(part1_(1)) + string(part1_(2));
            NASmcode = part1 + mcodeParts(2) + mcodeParts(3);
        else
            mcodeParts = split(mcode, "-");
            NASmcode = mcodeParts(2) + mcodeParts(3);
        end
        %         sourceDir = sprintf("/run/user/1664601583/gvfs/smb-share:server=badura_nas.local,share=data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
        sourceDir = sprintf("/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s", NASmcode);
        %     sourceDir = sprintf("/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
        %     sourceDir = sprintf("/run/user/1664601583/gvfs/smb-share:server=badura_nas.local,share=data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
        destDir = sprintf('/home/i.klinkhamer/Documents/Data/behaviorData/%s/', mcode);


        % List all files in the source directory
        folderList = dir(fullfile(sourceDir));
        if isempty(folderList)
            sourceDir = sprintf("/run/user/1664601583/gvfs/smb-share:server=badura_nas.local,share=data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
            folderList = dir(fullfile(sourceDir));
        end
        if isempty(folderList)
            sourceDir = sprintf("/run/user/1664601583/gvfs/smb-share:server=badura_nas.local,share=data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/%s/", NASmcode);
            folderList = dir(fullfile(sourceDir));
        end
        if isempty(folderList)
            disp("Folders for mouse not found. Check path to NAS.")
            disp("To find path to NAS, go to NAS and drag a file into matlab. Then copy the path from the load function in the terminal.")
        end

        folders = "";
        idx = 1;
        ses = 1;
        for i = 1 : length(folderList)
            if folderList(i).isdir && ~strcmp(folderList(i).name, '.') && ~strcmp(folderList(i).name, '..')
                sessionDates(m).sessions(ses).s = folderList(i).name;
                ses = ses + 1;
                folders(idx) = string(fullfile(folderList(i).folder, folderList(i).name));
                idx = idx + 1;
            end
        end

        for j = startFolder:length(folders)
            secondSessionOfTheDay = false;
            fileList = dir(folders(j));
            folderExtensions = "";
            for jj = 1 : length(fileList)
                [~, ~, folderExt] = fileparts(fileList(jj).name);
                folderExtensions(jj) = string(folderExt);
            end
            fileList = fileList(~[fileList.isdir] & folderExtensions == ".mat"); % You can specify the file extension you want to copy

            skipfile = zeros(1,length(fileList));
            lastSession = 1;
            fileSessions = [];
            for ii = 1:numel(fileList)
                parts = split(fileList(ii).name, '_');
                pattern = 's\d+\d+';
                mask = regexp(parts, pattern);
                matching_idx = find(cellfun(@(x) ~isempty(x), mask));
                if isempty(matching_idx)
                    skipfile(ii) = true;
                    continue
                end
                parts2 = split(parts(matching_idx), "");
                fileSessions(ii) = str2double(string(parts2{3}) + string(parts2{4}));
                lastSession = max(lastSession, str2double(string(parts2{3}) + string(parts2{4})));
            end
            fileSessions = unique(fileSessions);

            try
                partsFile = split([fileList.name], '_');
            catch
                disp("You need to connect the NAS")
            end
            pattern = 's\d+\d+';
            maskPattern = regexp(partsFile, pattern);
            matching_indices = find(cellfun(@(x) ~isempty(x), maskPattern));
            newMaskPatternFiles = partsFile(matching_indices);
            newMaskFilestest = [];
            for p = 1:length(fileSessions)
                pSes = fileSessions(p);
                if pSes < 10
                    patternNewtest = sprintf("s0%d", pSes);
                else
                    patternNewtest = sprintf("s%d", pSes);
                end
                newMaskFilestest(p).files = strcmp(newMaskPatternFiles, patternNewtest);
            end
            num_files_session = arrayfun(@(x) sum(x.files), newMaskFilestest);
            completed_sessions_mask = num_files_session >= 482 & num_files_session <= 483 | num_files_session >= 723;
            %         completed_sessions_mask = arrayfun(@(x) sum(x.files) >= 482 & sum(x.files) <= 483 | sum(x.files) >= 723, newMaskFilestest);
            if sum(completed_sessions_mask) > 1
                secondSessionOfTheDay = true;
                %             disp("test")
            end
            completedSessions = fileSessions(completed_sessions_mask);
            if ~isempty(completedSessions)
                if completedSessions(1) < 10
                    newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(1)));
                else
                    newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s%d", completedSessions(1)));
                end
                if completedSessions(end) < 10
                    newMaskFiles2 = strcmp(newMaskPatternFiles, sprintf("s0%d", completedSessions(end)));
                else
                    newMaskFiles2 = strcmp(newMaskPatternFiles, sprintf("s%d", completedSessions(end)));
                end
            else
                [~, idx] = max(num_files_session); % IK TODO: now it just takes the biggest session, but maybe add a way to combine multiple incompleted sessions.
                newMaskFiles = strcmp(newMaskPatternFiles, sprintf("s%d", fileSessions(idx)));
            end

            if length(newMaskFiles) < length(fileList)
                nExtraIdcs = length(fileList) - length(newMaskFiles);
                newMaskFiles(length(newMaskFiles)+1:length(newMaskFiles)+nExtraIdcs) = 0;
                newMaskFiles2(length(newMaskFiles2)+1:length(newMaskFiles2)+nExtraIdcs) = 0;
            end

            fprintf('|----'); % Print initial bar
            prevPercent = 0;

            % Loop through each file and copy it to the destination directory
            for i = 1:numel(fileList)
                behaviorDataResultsDir = "/home/i.klinkhamer/Documents/Data/behaviorDataResults/";
                if newMaskFiles(i)
                    [~, folderName, ~] = fileparts(folders(j));
                    destFolder = fullfile(destDir, folderName);
                    if ~exist(destFolder, 'dir')
                        mkdir(destFolder);
                    end

                    behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderName);
                    if ~exist(behaviorDataResultsFolder, 'dir')
                        mkdir(behaviorDataResultsFolder);
                    end
                    sourceFolder = fullfile(sourceDir, folderName);
                    sourceFile = fullfile(sourceFolder, fileList(i).name);
                    destFile = fullfile(destFolder, fileList(i).name);
                    if ~exist(destFile, "file")
                        copyfile(sourceFile, destFile);
                    end
                end
                if secondSessionOfTheDay && newMaskFiles2(i)
                    [~, folderName, ~] = fileparts(folders(j));
                    folderNameNew = sprintf("%s2", folderName);
                    destFolder = fullfile(destDir, folderNameNew);
                    if ~exist(destFolder, 'dir')
                        mkdir(destFolder);
                    end

                    behaviorDataResultsFolder = fullfile(behaviorDataResultsDir, mcode, folderNameNew);
                    if ~exist(behaviorDataResultsFolder, 'dir')
                        mkdir(behaviorDataResultsFolder);
                    end
                    sourceFolder = fullfile(sourceDir, folderName);
                    sourceFile = fullfile(sourceFolder, fileList(i).name);
                    destFile = fullfile(destFolder, fileList(i).name);
                    if ~exist(destFile, "file")
                        copyfile(sourceFile, destFile);
                    end
                end

                numValidFiles = numel(fileList(newMaskFiles | (secondSessionOfTheDay & newMaskFiles2)));
                percentComplete = i / numValidFiles * 100;

                if floor(percentComplete) > prevPercent || percentComplete == 100
                    fprintf(repmat('\b', 1, numel(num2str(floor(percentComplete))) + 3));
                    fprintf('*');
                    fprintf('| %d%%', floor(percentComplete));
                    prevPercent = floor(percentComplete) + 1;
                end

            end
            fprintf("folder %d done \n", j)
        end

        batchProcessTrialsExampleIK(mcode)
        try
            rmdir(destDir, 's'); % Only remove the DESTdir!!! If you remove the SOURCEdir you delete all the data from the NAS!!!
        catch
        end
    end

    fprintf("\nLoop done, waiting for new UI input confirmation.\n")
end
end
