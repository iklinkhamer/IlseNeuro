% Reexport simpleSpike rasters (via JK_PSTH_Selection) of the following
% units:

function tmp_psth_reexport_switchRasters_IK(mouseName, kwargs)
arguments
    mouseName = "Ana2";
    kwargs.saveFigs = true;
    kwargs.outputFolder = fullfile(Env.getBayesLabUserRoot,"/TraceExperiments/AnalysisOutput/Trace C4 Figures/switch_cells/rasters 10% contamination good units");
    kwargs.directory = fullfile(Env.getBayesLabUserRoot, "/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
    kwargs.evaluateC4Analysis = true;
end

import JkUtils.time.dateInt2datetime;
import JkUtils.flatmap;

stitchingName = "SwitchSessionStitching";

if not(isfolder(kwargs.outputFolder))
    mkdir(kwargs.outputFolder)
end
switchSession = Subject(mouseName).collectStitchedEphysSessions(stitchingName);
units = switchSession.collectKilosortUnits();

if kwargs.evaluateC4Analysis
    % Define the file path
    file_path = fullfile(kwargs.directory, "SwitchSessionStitching", ...
        "c4", "c4_results_fpfnThreshold_0.1_confidenceRatio_1.5", "cluster_predicted_cell_type.tsv");

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
        warning("File 'cluster_predicted_cell-type.tsv' not found.");
    end

    % Display extracted neuron numbers
    disp(neuron_ids);


    if kwargs.saveFigs
        for t = 1:size(cell_types,1)
            cell_type = cell_types(t);
            idcs_neurons_of_this_cell_type = classification_neurons==cell_type;
            neurons_of_this_cell_type = neuron_ids(idcs_neurons_of_this_cell_type);
            neuron_specs = ...
                { mouseName neurons_of_this_cell_type
                };

            outputFolder = fullfile(kwargs.outputFolder, cell_type);
            cellfun ...
                ( @(name, clusters) ...
                Subject(name) ...
                .collectStitchedEphysSessions(stitchingName) ...
                .JK_PSTH_Selection_StitchedSessionsWrapper ...
                ( neuronFilterFn=@(neuron) ismember(neuron.index, clusters) ...
                , outputFolder=outputFolder ...
                , batchMode = true ...
                ) ...
                , neuron_specs(:,1), neuron_specs(:,2) ...
                , UniformOutput=false ...
                );

            %IK.IK_PSTH_Selection(units, outputFolder=fullfile(kwargs.outputFolder, cell_type), selectBatchMode=true, selectArray=neurons_of_this_cell_type)

        end
    end
    neuronIDs = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), [units.id]);
    [~,neuronIDs_filtered_in] = intersect(neuronIDs,neuron_ids);
    mask = true(1,length(neuronIDs));
    mask(neuronIDs_filtered_in) = false;
    neuronIDs_filtered_out_units = neuronIDs(mask);

    neuron_specs = ...
        { mouseName neuronIDs_filtered_out_units
        };

    outputFolder = fullfile(fileparts(kwargs.outputFolder), "rasters 10% contamination filtered-out units");
    cellfun ...
        ( @(name, clusters) ...
        Subject(name) ...
        .collectStitchedEphysSessions(stitchingName) ...
        .JK_PSTH_Selection_StitchedSessionsWrapper ...
        ( neuronFilterFn=@(neuron) ismember(neuron.index, clusters) ...
        , outputFolder=outputFolder ...
        , batchMode = true ...
        ) ...
        , neuron_specs(:,1), neuron_specs(:,2) ...
        , UniformOutput=false ...
        );

    %IK.IK_PSTH_Selection(units, outputFolder=fullfile(fileparts(kwargs.outputFolder), "rasters 10% contamination filtered-out units"), selectBatchMode=true, selectArray=neuronIDs_filtered_out_units)


    

else

    % outputFolder = fullfile(fileparts(kwargs.outputFolder), "rasters (not filtered)");
    % Subject(name) ...
    %     .collectStitchedEphysSessions(stitchingName) ...
    %     .JK_PSTH_Selection_StitchedSessionsWrapper ...
    %     ( outputFolder=outputFolder ...
    %     , batchMode = false ...
    %     ) ...

    %IK.IK_PSTH_Selection(units, outputFolder=fullfile(fileparts(kwargs.outputFolder), "rasters (not filtered)"))

end

end


