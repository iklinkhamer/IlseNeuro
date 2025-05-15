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

function inspectMouseUnits_CSG(mouseName, kwargs)
arguments
    mouseName = "Georgetown_Cbx";
    kwargs.saveFigs = true;
    kwargs.evaluateC4Analysis = false;
    kwargs.batchMode=true;
    kwargs.outputFolder = fullfile(Env.getBayesLabUserRoot,"/ContextMouseExperiments/Ilse/CSG/AnalysisOutput/Rasters/");
    kwargs.directory = fullfile(Env.getBayesLabUserRoot, "/ContextMouseExperiments/Ilse/CSG/");   
    kwargs.data_directory = fullfile(Env.getBayesLabUserRoot, "/ContextMouseExperiments/Analysis/BEPAK/analysis_directories/")
    kwargs.batchSelection = 1;
end
getSessionsLater = 1;
%i= 0;

ephys_session_dir = dir(fullfile(kwargs.directory, sprintf('StructEphysData_%s*.mat', mouseName)));
ephys_session_files = string({ephys_session_dir.name});
% ephysDataVec_ = [];
% for f = 1:length(ephys_session_files)
%     csg_data = load(fullfile(kwargs.dataFolder, ephys_session_files(f)));
%     neuron_struct = csg_data.unit.neuron;
%     ephysDataVec_= [ephysDataVec_, neuron_struct];
% end

sessions = 1:length(ephys_session_files);

mdir = dir(fullfile(kwargs.data_directory, mouseName));

    folder_names = {mdir.name};
    numeric_suffix = cellfun(@(n) regexp(n, '\d+$', 'match', 'once'), folder_names, 'UniformOutput', false);
    numeric_vals = cellfun(@(s) str2double(s), numeric_suffix);
    if getSessionsLater
        sessions = numeric_vals(~isnan(numeric_vals));
    end
    valid_indices = [mdir.isdir] & ~ismember(folder_names, {'.', '..'}) & ismember(numeric_vals, sessions);
    mfolders = string({mdir(valid_indices).name});

for s=1:length(sessions)
    folder = mfolders(s);
    % session = sessions(:)'
   % session = sessions(s)';
    %i = i+1;
    %if i ~=2
    %     continue
    %end
    % [fig, ax] = viewRaster(neuron, "cs", trialSubset = 1:80, ax = rasterAxes(1)
    %   


    try % AnalyzedEphys might not exist
        csg_data = load(fullfile(kwargs.directory, ephys_session_files(s)));
        units = csg_data.unit;
        units.neuron = units.neuron(kwargs.batchSelection:end);
        % neurons = csg_data.unit.neuron;
        % stimtimes = csg_data.unit.stimtimes;
    catch       
        sprintf("Loading units failed for session %d", s)
        continue
    end

    [units.neuron(:).id] = deal(units.neuron(:).irc_trueID);


    for i = 1:numel(units.neuron)
        units.neuron(i).id = double(units.neuron(i).id);
    end
    [units.neuron(:).session_id] = deal(folder);


    if not(isfolder(kwargs.outputFolder))
        mkdir(kwargs.outputFolder)
    end

    if kwargs.evaluateC4Analysis

        % Define the file path
        file_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, kwargs.results_subfolder, "cluster_predicted_cell_type.tsv");

        % Initialize variables
        cell_types = [];
        classification_neurons = [];
        neuron_ids = [];

        % Check if the file exists
        if exist(file_path, 'file')
            % Read the TSV file
            data = readtable(file_path, 'FileType', 'text', 'Delimiter', '\t');

            % Validate required columns
            if all(ismember(["cluster_id", "predicted_cell_type"], data.Properties.VariableNames))
                % Extract neuron numbers
                neuron_ids = sort(unique(double(data.cluster_id)));

                % Extract cell types
                cell_types = unique(data.predicted_cell_type);

                % Store classification data
                classification_neurons = string(data.predicted_cell_type);

            else
                warning("Required columns 'cluster_id' and 'predicted_cell_type' not found in the file.");
            end
        else
            warning("File 'cluster_predicted_cell_type.tsv' not found.");
        end

        % Display extracted neuron numbers
        disp(neuron_ids);

        if kwargs.saveFigs
            for t = 1:size(cell_types,1)
                cell_type = cell_types(t);
                idcs_neurons_of_this_cell_type = classification_neurons==cell_type;
                neurons_of_this_cell_type = neuron_ids(idcs_neurons_of_this_cell_type);
                IK.IK_PSTH_Selection_CSG(mouseName, units, outputFolder=fullfile(kwargs.outputFolder, cell_type), subfolder = "rasters 10% contamination good units", selectBatchMode=true, selectArray=neurons_of_this_cell_type, folder = folder)

            end
        end
        neuronIDs = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), [units.id]);
        [~,neuronIDs_filtered_in] = intersect(neuronIDs,neuron_ids);
        mask = true(1,length(neuronIDs));
        mask(neuronIDs_filtered_in) = false;
        neuronIDs_filtered_out_units = neuronIDs(mask);

        IK.IK_PSTH_Selection_CSG(mouseName, units, outputFolder=kwargs.outputFolder, subfolder = "rasters 10% contamination filtered-out units", selectBatchMode=true, selectArray=neuronIDs_filtered_out_units, folder = folder)

    else
        if ~kwargs.batchMode
            IK.IK_PSTH_Selection_CSG(mouseName, units, outputFolder=kwargs.outputFolder, subfolder = "rasters (not filtered)", folder = folder)
        else
            neuronIDs = [units.neuron.id];
            % neuronIDs = neuronIDs_(selection);
            IK.IK_PSTH_Selection_CSG(mouseName, units, outputFolder=kwargs.outputFolder, subfolder = "rasters (not filtered)", selectBatchMode=true, selectArray=neuronIDs, folder = folder)
        end

    end
end
end


%%%
% classification_folder_path = fullfile(kwargs.directory, mouseName+"_"+session.timestampIdStr, kwargs.results_subfolder, "cell_type_classification");
%
%         neuron_numbers = [];    % Initialize an empty array to store numbers
%
%         files = dir(classification_folder_path);    % Get a list of all files in the target folder
%
%         for j = 1:length(files)
%             filename = files(j).name;
%
%             if startsWith(filename, "unit") % Check if the filename starts with "unit"
%                 neuron_number = regexp(filename, '\d+', 'match');   % Extract numbers from the filename using a regular expression
%                 if ~isempty(neuron_number)
%                     neuron_numbers = [neuron_numbers, str2double(neuron_number{1})]; %#ok<AGROW>
%                 end
%             end
%         end
%
%         neuron_numbers = sort(neuron_numbers);
%         %load(fullfile(classification_folder_path, "neurons_filtered_c4.mat"));
%         disp(neuron_numbers); % Outputs: [1 23 42]
% 
%         c4_folder = fileparts(classification_folder_path);
%         classification_neurons_ = readtable(fullfile(c4_folder,"cluster_predicted_cell_type.tsv"), 'FileType', 'text', 'Delimiter', '\t');
%         cell_types_ = unique(classification_neurons_(2:end,2));
%         cell_types = cell_types_.predicted_cell_type;
% 
%         classification_neurons = string(classification_neurons_.predicted_cell_type);
%         neuron_ids = double(classification_neurons_.cluster_id);
        %%%
