% Created on Mon Jan 13 11:43:52 2025
% Author: Ilse Klinkhamer

% Confusion Matrix Visualization Script
% 
% This script processes neural classification data from electrophysiology 
% experiments to compute and visualize the confusion matrix. It categorizes 
% neurons based on their baseline firing rates, compares true labels with 
% predicted labels from ensemble models, and plots both raw and normalized 
% confusion matrices for performance evaluation. 
%
% Dependencies:
% - Requires the npy-matlab package for loading .npy files.
% - Assumes custom classes for handling subjects, sessions, and neurons.

function confusionMatrix(mouseName, kwargs)
arguments
    mouseName = "Quimper"
    kwargs.neuronFilterFn = @(neuron) true;
    kwargs.outputFolder = "~/Documents/c4_neurons_temp_output/confusionMatrix";
    kwargs.subfolder54esdr56t2
    kwargs.fileTypes (1,:) string = ["png", "epsc"]; % "fig",
    kwargs.trialGroupingForPsth
    kwargs.printSize = [602 291]; % rectangular %[334 291] % square-ish aspect ratio for single psth plot
    kwargs.batchMode (1,1) logical = false % if true, save all units for which neuronFilterFn returns true
    kwargs.lineWidth (1,1) {isreal} = 0.5;
    kwargs.selectBatchMode (1,1) logical = false % if true, save a selection of neurons provided in an array
    kwargs.selectArray = [];
    kwargs.c4_folder = "c4";
    kwargs.c4_res_folder = "cell_type_classification";
    kwargs.labels = [1, 2, 3, 4, 5, 6; "GoC", "MLI", "MFB", "PkC_ss", "PkC_cs", "unlabeled"];
    kwargs.saveMismatchFigs = true;
end
close all

disp('Confusion Matrix Visualization Script');

mouse = Subject(mouseName);
sessions = mouse.collectSessions(sessionkinds="EPHYS");
labels = struct;
labels(1).names = kwargs.labels(2,:);
idcs = struct;


for s = 1:length(sessions)
    session = sessions(s);
    % raw_probabilities = readNPY(fullfile(session.ephysFolder, kwargs.c4_folder, "ensemble_predictions_ncells_nclasses_nmodels.npy")); % Requires npy-matlab package

    try
        neurons = session.collectKilosortUnits();
    catch
        disp("Analyzed ephys file probably not found")
        continue
    end

    baserates = neurons.baseRate();
    baselineFreqRanges = [0 3; 8 50; 50 inf; 3 8]; % Each row defines a range
    correspondingNeuronTypes = [5 2 4 6; "PkC_cs", "MLI", "PkC_ss", "unlabeled"];
    labels(s).frequency = zeros(size(baserates'));

    for i = 1:size(baselineFreqRanges, 1)
        type = str2double(correspondingNeuronTypes(1,i));
        labels(s).frequency(baserates >= baselineFreqRanges(i, 1) & baserates < baselineFreqRanges(i, 2)) = type;
    end
 
    % prediction_matrix = mean(raw_probabilities, 3); % Compute prediction matrix and predicted labels
    % [~, labels(s).predicted] = max(prediction_matrix, [], 2);

    cluster_predicted_cell_type = readmatrix(fullfile(session.ephysFolder, kwargs.c4_folder, "cluster_predicted_cell_type.tsv"), ...
        'Delimiter', '\t', 'OutputType', 'string', 'FileType', 'text');
    neuronIDs_c4 = str2double(cluster_predicted_cell_type(:,1));
    predicted_ = cluster_predicted_cell_type(:,2);
    for n = 1:length(predicted_)
        labels(s).predicted(n,1) = find(labels(1).names == predicted_(n));
    end
    neuronIDs_trace = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), [neurons.id]);
    [correspondingIDs, idcs(s).c4, idcs(s).trace] = intersect(neuronIDs_c4, neuronIDs_trace);
    freq_labels_sess = labels(s).frequency(idcs(s).trace);
    pred_labels_sess = labels(s).predicted(idcs(s).c4);

    wrongly_labeled_units_mask = freq_labels_sess ~= pred_labels_sess;
    wrongly_labeled_units = correspondingIDs(wrongly_labeled_units_mask);

    T = table(correspondingIDs, baserates(idcs(s).trace)',labels(1).names(freq_labels_sess)', labels(1).names(pred_labels_sess)', 'VariableNames', {'Neuron ID', 'baserate', 'frequency prediction', 'c4 prediction'});
    disp(T); % Display the Table
    writetable(T, fullfile(session.ephysFolder,kwargs.c4_folder, kwargs.c4_res_folder, 'confusionMatrixTable.csv'), 'Delimiter', '\t');

    if kwargs.saveMismatchFigs
        IK.IK_PSTH_Selection(neurons, outputFolder=fullfile(kwargs.outputFolder, "classifaction_mismatch_units"), selectBatchMode=true, selectArray=wrongly_labeled_units)
    end
    cm = confusionmat(freq_labels_sess, pred_labels_sess); % Compute confusion matrix

    disp('Confusion matrix, without normalization:'); % Display confusion matrix
    disp(cm);

    corresponding_labels = labels(1).names(unique([freq_labels_sess, pred_labels_sess]));

    if not(isfield(kwargs, "subfolder"))
        switch class(session)
            case 'Session'
                subfolder = session.tracePriorName();
            case 'StitchedSessions'
                subfolder = "";
            otherwise
                warning("Unkown session class %s", class(session))
                keyboard
        end
    else
        subfolder = kwargs.subfolder;
    end

    figure; % Plot confusion matrix
    plot_confusion_matrix(cm, corresponding_labels, 'Confusion Matrix', parula);
    fname = sprintf('confusionMatrix_%s', session.id);
    printFigure ...
        ( gcf ...
        , fname ...
        , folder = fullfile(kwargs.outputFolder, subfolder) ...
        , formats = kwargs.fileTypes ...
        , size = kwargs.printSize ...
        )

    cm_normalized = cm ./ sum(cm, 2); % Normalize confusion matrix
    disp('Normalized confusion matrix:');
    disp(cm_normalized);

    figure; % Plot normalized confusion matrix
    plot_confusion_matrix(cm_normalized, corresponding_labels, 'Normalized Confusion Matrix', parula);
    fname2 = sprintf('normalizedConfusionMatrix_%s', session.id);

    printFigure ...
        ( gcf ...
        , fname2 ...
        , folder = fullfile(kwargs.outputFolder, subfolder) ...
        , formats = kwargs.fileTypes ...
        , size = kwargs.printSize ...
        )
end
% Initialize empty arrays to store all labels
all_frequency_labels = [];
all_predicted_labels = [];

% Loop through each session
for i = 1:length(labels)
    % Extract frequency labels and predicted labels for the current session
    current_frequency_labels = labels(i).frequency(idcs(i).trace);
    current_predicted_labels = labels(i).predicted(idcs(i).c4);

    % Concatenate the labels to the overall arrays
    all_frequency_labels = [all_frequency_labels; current_frequency_labels];
    all_predicted_labels = [all_predicted_labels; current_predicted_labels];
end

% Compute the confusion matrix
cm = confusionmat(all_frequency_labels, all_predicted_labels);
corresponding_labels = labels(1).names(unique([all_frequency_labels, all_predicted_labels]));

% Display the confusion matrix
disp('Confusion matrix across all sessions:');
disp(cm);

figure; % Plot confusion matrix
plot_confusion_matrix(cm, corresponding_labels, 'Confusion Matrix', parula);
fname = sprintf('allSessionsConfusionMatrix_%s', mouseName);

printFigure ...
    ( gcf ...
    , fname ...
    , folder = kwargs.outputFolder ...
    , formats = kwargs.fileTypes ...
    , size = kwargs.printSize ...
    )


cm_normalized = cm ./ sum(cm, 2); % Normalize confusion matrix
disp('Normalized confusion matrix:');
disp(cm_normalized);

figure; % Plot normalized confusion matrix
plot_confusion_matrix(cm_normalized, corresponding_labels, 'Normalized Confusion Matrix', parula);
fname2 = sprintf('allSessionsNormalizedConfusionMatrix_%s', mouseName);

printFigure ...
    ( gcf ...
    , fname2 ...
    , folder = kwargs.outputFolder ...
    , formats = kwargs.fileTypes ...
    , size = kwargs.printSize ...
    )

end

function plot_confusion_matrix(cm, labels, title_text, cmap) % Function to plot confusion matrix
    imagesc(cm);
    colormap(cmap);
    colorbar;
    title(title_text, 'Interpreter', 'none');
    xlabel('Predicted Label');
    ylabel('Frequency Label');
    xticks(1:length(labels));
    yticks(1:length(labels));
    xticklabels(labels);
    yticklabels(labels);
    xtickangle(45);
    axis square;
    set(gca, 'TickDir', 'out');
end



    % 
    % baselineFreqRanges = [0 3; 3 8; 8 50; 50 inf]; 
    % correspondingNeuronTypes = ["PkC_cs", "unlabeled", "MLI", "PkC_ss"];
% 
% % Define 1 path
% dir_path = "/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/Seattle/Seattle_20200909133404/c4/cell_type_classification/";

% labels_ = [5, 4, 3, 2, 1; "PkC_cs", "PkC_ss", "MFB", "MLI", "GoC"]; 

% labels_ = readmatrix(fullfile(classificationFolder2, "label_correspondence.tsv"), ...
%                     'Delimiter', '\t', 'OutputType', 'string', 'FileType', 'text');
% labels_(:, end) = []; % Remove the last column

% import ComplexSpikeCuration.isPotentialComplexSpikeChannel
% isComplexNeuron = isPotentialComplexSpikeChannel(neurons(1), 5);


% directory_path = fullfile("/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
% directory = dir(directory_path);
% foldernames = {directory.name};
% ephysfolders_mask = [directory.isdir] & contains(foldernames, mouseName);
% ephysfolders_idcs = find(ephysfolders_mask);



% sessionFolder_idx = ephysfolders_idcs(s);
% sessionFolder = directory(sessionFolder_idx).name;
% classificationFolder = fullfile(directory_path, sessionFolder, classificationFolderName);
% classificationFolder2 = fullfile(directory_path, sessionFolder, classificationFolderName, "cell_type_classification");


