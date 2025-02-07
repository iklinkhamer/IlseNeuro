%% IK 26-4-24 Add NaN values to the trial data.
function behavior_trial_data = addNantoTrialData(mInput)
p = IkUtils.getParams();
if nargin < 1

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to process?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );
    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);

else
    if any(contains(mInput, p.mouseNames))
        mcode = p.mouseList(find(p.mouseNames == mInput));
    else
        mcode = mInput;
    end
end



path_list = getBehaviorPathList(mcode);


for f = 1 : length(path_list)

    folder = path_list(f);
    try
        load(fullfile(folder, 'trialdata.mat'));
    catch
        folder = path_list(f);
        s = "";
        filename = recursiveFileSearch("trialdata*.mat", mcode, s, folder);
        try
            load(fullfile(folder, filename(end)));
        catch
            continue
        end
    end
    continuefromnow = false;
    lastTrial = length(behavior_trial_data.eyelidpos);

    if lastTrial < 240
        trialdata2 = behavior_trial_data;
        for i = fieldnames(behavior_trial_data)'            
            var = trialdata2.(i{1});
            new_var = var;
            [rows, cols] = size(var); 
            if i == "meanAll"
                continuefromnow = true;
            end
            lastTrial = length(behavior_trial_data.eyelidpos);
            if i == "CRamp"
                continuefromnow = false;
                lastTrial = length(trialdata2.CRamp);
            end
            if continuefromnow
                continue
            end
            if isnumeric(var) && isequal(class(var), 'double') && i ~= "CRamp5"
                if rows == 1 || rows == 200
                    new_var = nan(rows, 240);
                    new_var(:,1:min(lastTrial, length(var))) = var;
                else
                    new_var = nan(240, cols);     
                    if i == "CRamp"
                        new_var = nan(220,cols);
                    end
                    new_var(1:min(lastTrial, size(var,1)),:) = var;
                end
%                 try
% 
%                     new_var(lastTrial,:) = nan;
%                 catch
%                     new_var(:,lastTrial) = nan;
%                 end
            elseif isequal(class(var), 'string') 
                new_var = var;
                new_var(lastTrial+1:240) = '';
            elseif isequal(class(var), 'cell')
                [rows, cols] = size(var);
%                 new_var = cell(rows, cols);
                %new_var{:,:} = nan;
                if rows == 1 || rows == 200   
                    new_var = cell(rows, 240);
                    for ii = 1:rows
                        for j = 1:lastTrial
                            new_var{ii, j} = var{ii,j};
                        end
                    end
                elseif cols == 1 || cols == 200 
                    new_var = cell(240, cols);
                    for ii = 1:min(lastTrial, length(var))
                        for j = 1:cols
                            new_var{ii, j} = var{ii, j};
                        end
                    end
                end

            end
            trialdata2.(i{1}) = new_var;
            %newfolder = fullfile(fileparts(folder),"230525");
            if i == "CRamp"
                continuefromnow = true;
            end
        end
        [~, folderName, ~] = fileparts(folder);
        if any(ismember(folderName, p.s))
            s = find(p.s == folderName);
        else
            s = folderName;
        end
        stamptrialor = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'originaltrialdata', file_extension = '.mat');
        save(fullfile(folder, stamptrialor), 'behavior_trial_data');
        behavior_trial_data_save = behavior_trial_data;
        behavior_trial_data = trialdata2;
        stamptrial = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialdata', file_extension = '.mat');
        save(fullfile(folder, stamptrial), 'behavior_trial_data');
        behavior_trial_data = behavior_trial_data_save;
        continuefromnow = false;
    end
end



end