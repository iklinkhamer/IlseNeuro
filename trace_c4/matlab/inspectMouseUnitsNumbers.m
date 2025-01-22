% INSPECTMOUSEUNITSNUMBERS Analyzes neuron units for a specified mouse.
% 
% This function inspects and processes neuron unit numbers for all
% electrophysiology (EPHYS) sessions of a given mouse. It collects unit
% information, outputs the count of neurons per session, and stores the
% results in a temporary output folder as a table with clear column separation.
% 
% USAGE:
%   inspectMouseUnitsNumbers(mouseName)
%
% INPUT:
%   mouseName (string, optional) - Name of the mouse to analyze. If not
%       provided, defaults to "Quimper".
%
% OUTPUT:
%   - Displays the number of sessions found and the number of neurons
%     detected per session.
%   - Outputs extracted neuron unit numbers for each session to the console.
%   - Saves the neuron counts per session in a table with tab-separated columns.
%
% NOTES:
%   - This script assumes specific folder structures for session files 
%     and neuron unit classification.
%   - Any missing "AnalyzedEphys" data or inaccessible directories will
%     trigger a warning and skip those sessions.
%
% DEPENDENCIES:
%   - Requires the Subject class and related methods like
%     `collectSessions` and `collectKilosortUnits`.
%   - Relies on `mkdir`, `dir`, and `writetable` for kwargs.directory and file management.
%
% Author: Ilse Klinkhamer
% Date: Mon Jan 13 11:43:52 2025

function inspectMouseUnitsNumbers(mouseName, kwargs)
arguments
    mouseName = "Iowa";
    kwargs.c4_folder = "c4";
    kwargs.outputFolder = "~/Documents/c4_neurons_temp_output/rasters 10% contamination good units";
    kwargs.directory = fullfile("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
end
mouse = Subject(mouseName);
sessions = mouse.collectSessions(sessionkinds="EPHYS");

fprintf("Found %d sessions for %s\n\n", numel(sessions), mouseName)

% Initialize variables
session_ids = {};
neuron_counts_c4 = [];
neuron_counts_trace = [];
neuron_counts_c4_old = [];

% Ensure the output folder exists
if not(isfolder(kwargs.outputFolder))
    mkdir(kwargs.outputFolder)
end

for session = sessions(:)'
    try
        % Load Kilosort units for the session
        units = session.collectKilosortUnits();
    catch
        fprintf("Warning: could not load AnalyzedEphys for session %s\n\n", session.timestampIdStr)
        continue
    end       
    session_id = session.timestampIdStr;  % Get the session ID
    % fprintf("Found %d neurons for session %s\n\n", numel(units), session_id);  
   
    try
        cluster_predicted_cell_type = readmatrix(fullfile(session.ephysFolder, kwargs.c4_folder, "cluster_predicted_cell_type.tsv"), ...
            'Delimiter', '\t', 'OutputType', 'string', 'FileType', 'text');
        neuronIDs_c4 = str2double(cluster_predicted_cell_type(:,1));        
    catch
        fprintf("Warning: classification folder does not exist for session %s\n", session_id);
    end

    numbers_c4_old = [];
    c_folder_paths = ["res_5perc_folder", "hide_res_folder", "res_5perc_old_figs_log_cwin2000/res_10perc_old_figs_log_cwin2000"];
    for f = 1:length(c_folder_paths)
        try
            classification_folder_path = fullfile(session.ephysFolder, "c4",c_folder_paths(f), "cell_type_classification");
            files = dir(classification_folder_path); % Get a list of all files in the target folder

            for j = 1:length(files)
                filename = files(j).name;

                % Check if the filename starts with "unit"
                if startsWith(filename, "unit")
                    % Extract numbers from the filename using a regular expression
                    number = regexp(filename, '\d+', 'match');
                    if ~isempty(number)
                        numbers_c4_old = [numbers_c4_old, str2double(number{1})]; %#ok<AGROW>
                    end
                end
            end
        catch
        end
    end
    
    session_ids{end+1} = session_id; %#ok<AGROW>   % Store session ID and neuron count
    neuron_counts_c4(end+1) = numel(neuronIDs_c4); %#ok<AGROW>
    neuron_counts_trace(end+1) = numel(units); %#ok<AGROW>
    neuron_counts_c4_old(end+1) = numel(numbers_c4_old); %#ok<AGROW>
    
    % disp(numbers_c4); % Display neuron numbers    % Outputs: [1 23 42]
end
results_table = table(session_ids', neuron_counts_c4_old',neuron_counts_c4',neuron_counts_trace', ...
    'VariableNames', {'SessionID', 'Neuron Count C4 Old Thres', 'Neuron Count C4','Neuron Count Trace'}); % Create a table with session IDs and neuron counts

output_file = fullfile(kwargs.outputFolder, mouseName + "_neuron_counts_per_session.tsv"); % Save the table to the output folder with tab-separated columns
writetable(results_table, output_file, 'FileType', 'text', 'Delimiter', '\t');

fprintf("Neuron counts table saved to: %s\n", output_file);


session_ids = cellstr(string(results_table.SessionID)); % Ensure session_ids is in the correct format (cell array of strings) % Convert to cell array of strings
neuron_counts = table2array(results_table(:, 2:end)); % Extract neuron count columns from the table % Extract numeric data (columns 2 to 4)

figure;% Create a grouped bar plot
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);% Set the figure to fullscreen
bar(categorical(session_ids), neuron_counts, 'grouped');

xlabel('Session ID', 'Interpreter', 'none'); % Add labels, title, and legend % Prevent LaTeX interpretation of underscores
ylabel('Neuron Counts');
title(['Neuron Counts Across All Categories for Mouse ' char(mouseName)]);
legend(results_table.Properties.VariableNames(2:end), 'Location', 'bestoutside');
xtickangle(45); % Rotate x-axis labels for better visibility
set(gca, 'TickDir', 'out');
grid on;

plot_output_file = fullfile(kwargs.outputFolder, mouseName + "_neuron_counts_comparison_plot.png"); % Save the plot as an image file
saveas(gcf, plot_output_file);

fprintf("Neuron counts comparison plot saved to: %s\n", plot_output_file);
end

% function inspectMouseUnitsNumbers(mouseName, kwargs)
% arguments
%     mouseName = "Quimper";
%     kwargs.c4_folder = "c4";
%     kwargs.outputFolder = "~/Documents/c4_neurons_temp_output";
% end
% kwargs.directory = fullfile("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
% 
% mouse = Subject(mouseName);
% sessions = mouse.collectSessions(sessionkinds="EPHYS");
% 
% fprintf("Found %d sessions for %s\n\n", numel(sessions), mouseName)
% 
% % Initialize variables
% session_ids = {};
% neuron_counts = [];
% 
% % Ensure the output folder exists
% if not(isfolder(kwargs.outputFolder))
%     mkdir(kwargs.outputFolder)
% end
% 
% for session = sessions(:)'
%     try
%         % Load Kilosort units for the session
%         units = session.collectKilosortUnits();
%     catch
%         fprintf("Warning: could not load AnalyzedEphys for session %s\n\n", session.timestampIdStr)
%         continue
%     end
% 
%     cluster_predicted_cell_type = readmatrix(fullfile(session.ephysFolder, kwargs.c4_folder, "cluster_predicted_cell_type.tsv"), ...
%         'Delimiter', '\t', 'OutputType', 'string', 'FileType', 'text');
%     neuronIDs_c4 = str2double(cluster_predicted_cell_type(:,1));
% 
%     % Get the session ID
%     session_id = session.timestampIdStr;
%     fprintf("Found %d neurons for session %s\n\n", numel(units), session_id);
% 
%     % Construct classification folder path
%     classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "cell_type_classification");
% 
%     % Initialize an empty array to store neuron numbers
%     numbers = [];
% 
%     % Get a list of all files in the target folder
%     if isfolder(classification_folder_path)
%         files = dir(classification_folder_path);
% 
%         for j = 1:length(files)
%             filename = files(j).name;
% 
%             % Check if the filename starts with "unit"
%             if startsWith(filename, "unit")
%                 % Extract numbers from the filename using a regular expression
%                 number = regexp(filename, '\d+', 'match');
%                 if ~isempty(number)
%                     numbers = [numbers, str2double(number{1})]; %#ok<AGROW>
%                 end
%             end
%         end
%     else
%         fprintf("Warning: classification folder does not exist for session %s\n", session_id);
%     end
% 
%     % Store session ID and neuron count
%     session_ids{end+1} = session_id; %#ok<AGROW>
%     neuron_counts(end+1) = numel(numbers); %#ok<AGROW>
% 
%     % Display neuron numbers
%     disp(numbers); % Outputs: [1 23 42]
% end
% 
% % Create a table with session IDs and neuron counts
% results_table = table(session_ids', neuron_counts', ...
%     'VariableNames', {'SessionID', 'NeuronCount'});
% 
% % Save the table to the output folder with tab-separated columns
% output_file = fullfile(kwargs.outputFolder, mouseName + "_neuron_counts_per_session.tsv");
% writetable(results_table, output_file, 'FileType', 'text', 'Delimiter', '\t');
% 
% fprintf("Neuron counts table saved to: %s\n", output_file);
% end



% % INSPECTMOUSEUNITSNUMBERS Analyzes neuron units for a specified mouse.
% % 
% % This function inspects and processes neuron unit numbers for all
% % electrophysiology (EPHYS) sessions of a given mouse. It collects unit
% % information, outputs the count of neurons per session, and stores the
% % results in a temporary output folder.
% % 
% % USAGE:
% %   inspectMouseUnitsNumbers(mouseName)
% %
% % INPUT:
% %   mouseName (string, optional) - Name of the mouse to analyze. If not
% %       provided, defaults to "Quimper".
% %
% % OUTPUT:
% %   - Displays the number of sessions found and the number of neurons
% %     detected per session.
% %   - Outputs extracted neuron unit numbers for each session to the console.
% %
% % NOTES:
% %   - This script assumes specific folder structures for session files 
% %     and neuron unit classification.
% %   - Any missing "AnalyzedEphys" data or inaccessible directories will
% %     trigger a warning and skip those sessions.
% %
% % DEPENDENCIES:
% %   - Requires the Subject class and related methods like
% %     `collectSessions` and `collectKilosortUnits`.
% %   - Relies on `mkdir` and `dir` for kwargs.directory management.
% %
% % Author: Ilse Klinkhamer
% % Date: Mon Jan 13 11:43:52 2025
% 
% function inspectMouseUnitsNumbers(mouseName)
% if nargin < 1 || isempty(mouseName)
%     mouseName = "Quimper"; % Default value
% end
% kwargs.outputFolder = "~/Documents/c4_neurons_temp_output";
% 
% kwargs.directory = fullfile("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
% 
% 
% 
% mouse = Subject(mouseName);
% 
% sessions = mouse.collectSessions(sessionkinds="EPHYS");
% 
% fprintf("Found %d sessions for %s\n\n", numel(sessions), mouseName)
% % i= 0;
% numbers_all_sessions = [];
% i = 0;
% 
% for session = sessions(:)'
%     i = i+1;
%     % i = i+1;
%     % if i ~=6
%     %     continue
%     % end
%     % 
% 
%    % try
% 
%         try % AnalyzedEphys might not exist
%             units = session.collectKilosortUnits();
%         catch
%             fprintf("Warning: could not load AnalyzedEphys for session %s\n\n", session.timestampIdStr)
%             continue
%         end
% 
%         fprintf("Found %d neurons for session %s\n\n", numel(units), session.timestampIdStr)
% 
%         if not(isfolder(kwargs.outputFolder))
%             mkdir(kwargs.outputFolder)
%         end
%         classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "cell_type_classification");
%         % classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "res_5perc_folder", "cell_type_classification");
%         % classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "hide_res_folder", "cell_type_classification");
%         % classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "res_5perc_old_figs_log_cwin2000", "res_10perc_old_figs_log_cwin2000", "cell_type_classification");
%         % Initialize an empty array to store numbers
%         numbers = [];
% 
%         % Get a list of all files in the target folder
%         files = dir(classification_folder_path);
% 
%         for j = 1:length(files)
%             filename = files(j).name;
% 
%             % Check if the filename starts with "unit"
%             if startsWith(filename, "unit")
%                 % Extract numbers from the filename using a regular expression
%                 number = regexp(filename, '\d+', 'match');
%                 if ~isempty(number)
%                     numbers = [numbers, str2double(number{1})]; %#ok<AGROW>
%                 end
%             end
%         end
% 
%         % Sort the numbers
%         numbers = sort(numbers);
%         numbers_all_sessions(i) = numel(numbers);
%         %load(fullfile(classification_folder_path, "neurons_filtered_c4.mat"));
%         disp(numbers); % Outputs: [1 23 42]
% 
% 
%         %IK.IK_PSTH_Selection(units, kwargs.outputFolder=kwargs.outputFolder, selectBatchMode=true, selectArray=numbers)
%    % catch
%    % end
% end
% numbers_all_sessions
% end

