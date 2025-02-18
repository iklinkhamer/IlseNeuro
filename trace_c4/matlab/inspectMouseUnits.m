% INSPECTMOUSEUNITS Analyzes neuron units for a given mouse across EPHYS sessions.
%
% This function processes electrophysiology (EPHYS) sessions for a specified
% mouse, identifying and analyzing neuron units. The function attempts to 
% extract unit information from session data and display the unit numbers. 
% Results are saved to a temporary output folder, and missing or inaccessible 
% session data is skipped with a warning.
%
% USAGE:
%   inspectMouseUnits(mouseName)
%
% INPUT:
%   mouseName (string, optional) - Name of the mouse to analyze. If not
%       provided, defaults to "Seattle".
%
% FUNCTIONALITY:
%   - Iterates through EPHYS sessions of the specified mouse.
%   - Collects neuron unit numbers and displays them for each session.
%   - Saves results to the specified output folder.
%   - Skips sessions with missing "AnalyzedEphys" data or inaccessible directories.
%
% DEPENDENCIES:
%   - Requires the Subject class with `collectSessions` and `collectKilosortUnits`.
%   - Relies on the `IK.IK_PSTH_Selection` function for processing units.
%
% OUTPUT:
%   - Displays the number of neurons for each session and their unit IDs.
%   - Processes and saves relevant data to the output folder.
%
% NOTES:
%   - Ensure proper folder structure for session and classification files.
%   - The function handles missing data gracefully but may skip certain sessions.
%
% Author: Ilse Klinkhamer
% Date: Mon Jan 13 11:43:52 2025

function inspectMouseUnits(mouseName, kwargs)
arguments
    mouseName = "Venice";
    kwargs.saveFigs = true;
    kwargs.outputFolder = fullfile(Env.getBayesLabUserRoot,"/TraceExperiments/AnalysisOutput/Trace C4 Figures/rasters 10% contamination good units");
    kwargs.directory = fullfile(Env.getBayesLabUserRoot, "/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
    kwargs.evaluateC4Analysis = false;
end

mouse = Subject(mouseName);

sessions = mouse.collectSessions(sessionkinds="EPHYS");

fprintf("Found %d sessions for %s\n\n", numel(sessions), mouseName)
%i= 0;

for session = sessions(:)'
    %i = i+1;
    %if i ~=2
    %     continue
    %end
    %



    try % AnalyzedEphys might not exist
        units = session.collectKilosortUnits();
    catch
        fprintf("Warning: could not load AnalyzedEphys for session %s\n\n", session.timestampIdStr)
        continue
    end

    fprintf("Found %d neurons for session %s\n\n", numel(units), session.timestampIdStr)

    if not(isfolder(kwargs.outputFolder))
        mkdir(kwargs.outputFolder)
    end

    if kwargs.evaluateC4Analysis
        classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, "c4", "cell_type_classification");

        neuron_numbers = [];    % Initialize an empty array to store numbers

        files = dir(classification_folder_path);    % Get a list of all files in the target folder

        for j = 1:length(files)
            filename = files(j).name;

            if startsWith(filename, "unit") % Check if the filename starts with "unit"
                neuron_number = regexp(filename, '\d+', 'match');   % Extract numbers from the filename using a regular expression
                if ~isempty(neuron_number)
                    neuron_numbers = [neuron_numbers, str2double(neuron_number{1})]; %#ok<AGROW>
                end
            end
        end

        neuron_numbers = sort(neuron_numbers);
        %load(fullfile(classification_folder_path, "neurons_filtered_c4.mat"));
        disp(neuron_numbers); % Outputs: [1 23 42]

        c4_folder = fileparts(classification_folder_path);
        classification_neurons_ = readtable(fullfile(c4_folder,"cluster_predicted_cell_type.tsv"), 'FileType', 'text', 'Delimiter', '\t');
        cell_types_ = unique(classification_neurons_(2:end,2));
        cell_types = cell_types_.predicted_cell_type;

        classification_neurons = string(classification_neurons_.predicted_cell_type);
        neuron_ids = double(classification_neurons_.cluster_id);



        if kwargs.saveFigs
            for t = 1:size(cell_types,1)
                cell_type = cell_types(t);
                idcs_neurons_of_this_cell_type = classification_neurons==cell_type;
                neurons_of_this_cell_type = neuron_ids(idcs_neurons_of_this_cell_type);
                IK.IK_PSTH_Selection(units, outputFolder=fullfile(kwargs.outputFolder, cell_type), selectBatchMode=true, selectArray=neurons_of_this_cell_type)

            end
        end
        neuronIDs = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), [units.id]);
        [~,neuronIDs_filtered_in] = intersect(neuronIDs,neuron_numbers);
        mask = true(1,length(neuronIDs));
        mask(neuronIDs_filtered_in) = false;
        neuronIDs_filtered_out_units = neuronIDs(mask);

        IK.IK_PSTH_Selection(units, outputFolder=fullfile(fileparts(kwargs.outputFolder), "rasters 10% contamination filtered-out units"), selectBatchMode=true, selectArray=neuronIDs_filtered_out_units)

    else
        IK.IK_PSTH_Selection(units, outputFolder=fullfile(fileparts(kwargs.outputFolder), "rasters (not filtered)"))

    end
end
end

