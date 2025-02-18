%% I.K. 17-02-2025 initiate spike sorting with ironclust
function initSpikeSorting2()

disp("Calculating...")

mousecodes = "Venice";
sessions = 1;

for mcode = mousecodes
    for s = sessions
        fprintf("Mouse: %s Session: %d\n", mcode, s)
        
        % dataFileTypes = ["RAWdata*.bin", "continuous*.dat", "Data4KS2.bin"];
        % for d = dataFileTypes
        %     [files, folders] = recursiveFileSearch(d, mcode, s);
        %     if ~isempty(files)
        %         dataPath = folders(1);
        %         fileName = files(1);
        %         break
        %     end
        % end
        % if ~exist("dataPath", "var")
        %     disp("data not found")
        %     if IkUtils.do_prompt_yes_or_no("Convert raw ephys data?")
        %         convertData()
        %     end
        %     return
        % end
        % 
        % config_file_pattern = "Config_h32_oe*.prm";
        % [files, folders] = recursiveFileSearch(config_file_pattern, mcode, s);
        % if isempty(files)
        %     stamp = nameDateStampFiles(mcode = mcode, s = s);
        %     copyfile("/home/i.klinkhamer/Documents/Config_h32_oe.prm", fullfile(P.pathSpikeSortingHome, mcode, P.s(s), sprintf("Config_h32_oe_%s.prm", stamp)));
        %     kwargs.setDataPath = true;
        %     [files, folders] = recursiveFileSearch(config_file_pattern, mcode, s);
        % end
        % for c = files
        %     configFileName_ = fullfile(mcode, P.s(s), c);
        %     if isfile(fullfile(which(configFileName_)))
        %         configFileName = fullfile(which(configFileName_));
        %         break
        %     end
        % end
        % if ~exist("configFileName", "var")
        %     disp("Params file not found")
        %     return
        % end
        % 
        % if resetProbePath
        %     probePath = which("DBC_3.1-64-H2_IK.prb");
        %     if isempty(probePath)
        %         probePath = P.prbPath;
        %     end
        %     setProbePath(configFileName, probePath);
        % end

        configFileName = "'/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Venice/Venice_20240523111730(copy)/Extraction2Bin/Config_h32_oe_2.prm'";
        configFileName = "'/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Venice/Venice_20240523111730(copy)/c4/continuous/Data_AP_LFP/Config_h32_oe_2.prm'";
        if isfile(erase(configFileName, ".prm") + "_full.prm")
            str = "irc_IK manual " + configFileName;
        else
            str = "irc_IK full " + configFileName;
        end

        eval(str)            
    end
end

disp("End of loop. Waiting for UI input.")

end

function setParamDataPath(filename, dataPath, fileName)
prm_file  = fopen(string(filename),'r');
prm_file_contents = fread(prm_file, 'char*1');
prm_file_contents_str = char(prm_file_contents)'; % Convert the file contents to a string
raw_data_var = 'vcFile = ';
probe_file_idx = strfind(prm_file_contents_str, raw_data_var); % Find the position of "probe_file" in the file contents
semicolon_idx = strfind(prm_file_contents_str(probe_file_idx:end), ';');  % Find the position of the semicolon after "probe_file"
contents_before = prm_file_contents_str(probe_file_idx + length(raw_data_var):probe_file_idx + semicolon_idx - 2); % Extract the substring between "probe_file" and the semicolon

prm_file_contents_new = strrep(prm_file_contents, sprintf("%s%s", raw_data_var, contents_before), sprintf("%s'%s'", raw_data_var, fullfile(dataPath, fileName)));

prm_file = fopen(filename,'w');
fprintf(prm_file, '%s', prm_file_contents_new);
fclose(prm_file);
end

function setProbePath(filename, probePath)
prb_file  = fopen(filename,'r');
prm_file_contents = fread(prb_file, 'char*1');
prm_file_contents_str = char(prm_file_contents)'; % Convert the file contents to a string
prb_file_var = 'probe_file = ';
probe_file_idx = strfind(prm_file_contents_str, prb_file_var); % Find the position of "probe_file" in the file contents
semicolon_idx = strfind(prm_file_contents_str(probe_file_idx:end), ';');  % Find the position of the semicolon after "probe_file"
contents_before = prm_file_contents_str(probe_file_idx + length(prb_file_var):probe_file_idx + semicolon_idx - 2); % Extract the substring between "probe_file" and the semicolon

prm_file_contents_new = strrep(prm_file_contents, sprintf("%s%s", prb_file_var, contents_before), sprintf("%s'%s'", prb_file_var, probePath));
prb_file = fopen(filename,'w');
fprintf(prb_file, '%s', prm_file_contents_new);
fclose(prb_file);
end

