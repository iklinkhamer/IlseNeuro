% Select (and export) rasters and PSTHs from the supplied set of neurons
%
%   trialGroupingForPsth :: a cell array specifying which trial labels to be plotted
%   together
%       Example: If the possible trial labels are [0 1 2 3 4 5 6 7] and the
%       trialGrouping is {[0 1] [2 3 4 5 6 7]}, then the PSTHs for trials labeled 0
%       and 1 are plotted together and the PSTHs for the trials labeled 2 through 7
%       are plotted together in different subplots.
%
%
% Copyright © 2024 NarainLab
% Erasmus Medical Center, Rotterdam, The Netherlands
% v2024-01-30
%
function IK_PSTH_Selection(neurons, kwargs)
arguments
    neurons (1,:) KilosortUnit
    % kwargs.neuronFilterFn = mkBaselineFilter(10,50)%mkMinBaselineFilter(50,inf)
    kwargs.neuronFilterFn = @(neuron) true;
    kwargs.outputFolder = fullfile( Env.getBayesLabUserRoot ...
        , "results" ...
        , "manuscripts" ...
        , "TRACE Complex Spike" ...
        , "matlabVisuals" ...
        , "simpleSpikeRasters" ...
        );
    kwargs.subfolder
    kwargs.fileTypes (1,:) string = ["png", "epsc"]; % "fig",
    kwargs.trialGroupingForPsth
    kwargs.printSize = [602 291]; % rectangular %[334 291] % square-ish aspect ratio for single psth plot
    kwargs.batchMode (1,1) logical = false % if true, save all units for which neuronFilterFn returns true
    kwargs.lineWidth (1,1) {isreal} = 0.5;
    kwargs.selectBatchMode (1,1) logical = false % if true, save a selection of neurons provided in an array
    kwargs.selectArray = [];
end

DnDisp('----------- Single Unit Analysis ----------')

if isfield(kwargs, 'trialGroupingForPsth')
    nPsthPlots = numel(kwargs.trialGroupingForPsth);
else
    nPsthPlots = 1;
end

nRasterPlots = 1;
nPlotsPerAlignment = nRasterPlots+nPsthPlots;
nAlignments = numel(["CS", "US"]);


[axs, fig] = JkUtils.initPlots([nPlotsPerAlignment nAlignments]);
rasterAxes = axs(1,:);
psthAxes = axs(2:end,:);

%% FUTURE paper -> square aspect ratio

% printSize = [334         291];

%% Regular
singleScreenPos = [86.3333  495.6667  602.0000  290.6667];
multiScreenPos = [1685 1138 602 291]; % regular
% printSize = [602 291]; -> kwargs.printSize
%% sessionStitching
% printSize = [398         447];
% multiScreenPos = [ 1         784        1049        1043];
% singleScreenPos = [ 1         784        1049        1043];

% fig.Position = [ 1         784        1049        1043];
switch JkUtils.hostname
    case {'bumblebee'}
        fig.Position = singleScreenPos;
    otherwise
        fig.Position = multiScreenPos;
end

nNeurons = numel(neurons)

for n = 1:nNeurons

    neuron = neurons(n);

    if not(kwargs.neuronFilterFn(neuron))
        % If a specific filtering function has been supplied (e.g. to check
        % for firing rate baseline) use it to skip neurons failing its
        % constraint.
        fprintf("\tSkipping neuron: %s\n", neuron.id)
        continue
    end

    neuron.viewRaster("cs", ax = rasterAxes(1))
    neuron.viewRaster ...
        ( "us" ...
        , ax = rasterAxes(2) ...
        , showTitle=false ...
        , showAxLabels=false ...
        )

    if isfield(kwargs, "trialGroupingForPsth")
        groupingArg = {"trialGroupingForPsth", kwargs.trialGroupingForPsth};
    else
        groupingArg = {};
    end
    neuron.viewPsthPerTrialType ...
        ( "cs" ...
        , groupingArg{:} ...
        , axs = psthAxes(:,1) ...
        , showTitle=false ...
        , showAxLabels=true ...
        , lineWidth = kwargs.lineWidth ...
        );
    neuron.viewPsthPerTrialType ...
        ( "us" ...
        , groupingArg{:} ...
        , axs = psthAxes(:,2) ...
        , showTitle=false ...
        , showAxLabels=true ...
        , lineWidth = kwargs.lineWidth ...
        );

    arrayfun(@hideAxesLabels, [psthAxes(1:end-1, 1); psthAxes(:,2)])

    % keyboard
    % linkaxes(psthAxes, 'y')
    syncY(psthAxes(:));

    % pause;
    % Find the last underscore
    neuronId = string(neuron.id);
    % Find the index of the last underscore
    lastUnderscoreIndex = strfind(neuronId, '_'); % Find all underscores
    lastUnderscoreIndex = lastUnderscoreIndex(end); % Get the last underscore
    % Extract the substring after the last underscore
    lastNumberStr = extractAfter(neuronId, lastUnderscoreIndex);
    % Convert the extracted substring to a number
    lastNumber = str2double(lastNumberStr);

    if ~kwargs.selectBatchMode
        saveFig = str2double(input('Want this one?','s'));
    else
        saveFig = 0;
    end
    if kwargs.batchMode || (kwargs.selectBatchMode && ismember(lastNumber, kwargs.selectArray)) || saveFig
        if not(isfield(kwargs, "subfolder"))
            switch class(neuron.session)
                case 'Session'
                    subfolder = neuron.session.tracePriorName();
                case 'StitchedSessions'
                    subfolder = "";
                otherwise
                    warning("Unkown session class %s", class(neuron.session))
                    keyboard
            end
        else
            subfolder = kwargs.subfolder;
        end


        fprintf("\tSaving neuron:%s\n", neuron.id)
        fname = sprintf('PSTHPair_%s', neuron.id);
        % filepath = fullfile( kwargs.outputFolder ...
        %                    , subfolder ...
        %                    , fname );
        printFigure ...
            ( fig ...
            , fname ...
            , folder = fullfile(kwargs.outputFolder, subfolder) ...
            , formats = kwargs.fileTypes ...
            , size = kwargs.printSize ...
            )
        % for filetype = kwargs.fileTypes(:)'
        % saveas(fig,filepath,filetype);
        % end
    end

    arrayfun(@cla, axs)
    arrayfun(@(ax) set(ax, 'YLimMode', 'auto'), axs)

end

% end

JkUtils.safeClose(fig)

end

function fn = mkBaselineFilter(mn, mx)

fn = @baselineFilter;

    function keep = baselineFilter(neuron)
        baseline = neuron.baseRate();
        keep = mn < baseline && baseline <= mx;
    end

end

function axs = syncY(axs)

lims = [axs.YLim];
lowers = lims(1:2:end-1);
uppers = lims(2:2:end);

newLims = [min(lowers) max(uppers)];

arrayfun(@(ax) ylim(ax, newLims), axs)

yTicks = axs(1).YTick;
set(axs, "YTick", yTicks)
set(axs, "YTickMode", "auto")
end

function hideAxesLabels(ax)
ax.XLabel.Visible = false;
ax.YLabel.Visible = false;

end
