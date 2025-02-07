%% IK & Nynke 2-11-23
function batchProcessTrialsExample(mcode)
addpath("/home/mick/Desktop/Ilse/ephys_code/");
addpath("/home/mick/Desktop/Ilse/ephys_code/Complex spike suite revised/");
addpath('/media/mick/DATA/Ilse/Data/convertedData/')
addpath("/home/mick/Desktop/Ilse/ephys_code/Code Nynke/analysis code/")
addpath("/home/mick/Desktop/Ilse/ephys_code/Code Nynke/neuroblinks-master/utilities/")
%base_dir = 'C:\Users\nmhet\Documents\studie\MEP\data';
if nargin < 1
%     mcode = "20-MI19442-03";

fprintf("\n-------------------------------------")
fprintf("\nWhich mouse would you like to curate?")
fprintf("\n-------------------------------------\n")

mouseCodes = arrayfun ...          % IK change
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
mouseNames = arrayfun ...
    ( @(mouse) string(mouse.name) ...
    , defaultMice() ...
    );
%mouseNames = [mouseNames, "All of the above"];


[mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
mcode = mouseCodes(mname_idx);
%mname = mouseNames(mname_idx);

end
p = getParams();


dir_content = dir(fullfile(p.pathBehaviorDataDATA, mcode));
path_list = IkUtils.getPathFromDir(dir_content);
% dates = dir_content(3:end).name;

for f = 1 : length(path_list)

    folder = path_list(f);

    if ~exist(fullfile(folder, 'trialdata.mat'), 'file')
        behavior_trial_data = processTrials(folder, 'recalibrate'); %,'recalibrate');  % Recalibrate eyelid
    else
        continue
    end


    if ~isempty(behavior_trial_data)

        save(fullfile(folder, 'trialdata.mat'), 'behavior_trial_data');

        baseline = mean(mean(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,1:40)));
        fullclosure = mean(max(behavior_trial_data.eyelidpos(behavior_trial_data.c_csdur==0,40:95),[],2));

        traces = (behavior_trial_data.eyelidpos - baseline) / (fullclosure - baseline);

        figure
        plot(behavior_trial_data.tm(1,:), traces)
    end
end
