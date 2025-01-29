%% I.K. 1-6-24
function channelPieChart(kwargs)
arguments
    kwargs.mouseTypesDefault (1,:) string = ["Shank2KOMUT", "Shank2KOWT"];
    kwargs.saveFigures (1,1) logical = true;
end
close all
app = getUImain(scriptname=mfilename);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)        
        break; % If the UI figure has been closed, exit the loop
    end
    disp("Calculating...")
    appData = evalin('base', 'appData');

    mouseCodes = appData.mousecodes;
    mouseTypes = appData.mousetypes;
    gene = appData.gene; loc = appData.loc; type = appData.type;
    if appData.loc == "N/A"
        Types = [gene+"KO"+"WT", gene+"KO"+"MUT", gene+"L7"+"KO", gene+"L7"+"MUT"];
    else
        Types = [gene+loc+"WT", gene+loc+"MUT"];
    end

    WTmice = mouseCodes(contains(mouseTypes, gene+loc+"WT"));
    MUTmice = mouseCodes(contains(mouseTypes, gene+loc+"MUT"));

    mice = {WTmice, MUTmice};

    P = IkUtils.getParams();
    outputRoot = fullfile ...
        (fullfile(P.figPath, 'pieCharts'));

    if ~exist(outputRoot, 'dir')
        mkdir(outputRoot)
    end

    for m = 1:length(Types)
        mouseCodesGroup = mice{m};

        Fig = pieChartForType(Types(m), false, mouseCodesGroup);
        Fig2 = pieChartForType(Types(m), false, mouseCodesGroup, 'cs', true);
        Fig3 = pieChartForType(Types(m), false, mouseCodesGroup, 'us', true);

        if kwargs.saveFigures
            saveas( Fig...
                , fullfile(outputRoot, sprintf("pieChart_%s.png", Types(m)))...
                );

            folder = outputRoot;

            file = sprintf("pieChart_%s.eps", Types(m));
            print(Fig, '-depsc', '-painters', fullfile(folder, file))

            saveas( Fig2...
                , fullfile(outputRoot, sprintf("pieChart_CS_%s.png", Types(m)))...
                );
            folder = outputRoot;
            file = sprintf("pieChart_CS_%s.eps", Types(m));
            print(Fig2, '-depsc', '-painters', fullfile(folder, file))

            saveas( Fig3...
                , fullfile(outputRoot, sprintf("pieChart_US_%s.png", Types(m)))...
                );
            folder = outputRoot;
            file = sprintf("pieChart_US_%s.eps", Types(m));
            print(Fig3, '-depsc', '-painters', fullfile(folder, file))
        end
    end
    disp("Loop end. Waiting for new UI input...")
end
end
function [fig] = pieChartForType(type, withExplode, mouseCodes, stim, general)

if nargin < 2
    withExplode = false;
end
if nargin < 4
    general = false;
end

fig = figure;

if general == false
    [ none ...
        , cs_facilitation, us_facilitation ...
        , cs_suppression, us_suppression ...
        , csus_facilitation, csus_suppression ...
        , cs_fac_us_sup ...
        , us_fac_cs_sup ] = countChannelsForType(type, mouseCodes);

    pieData = ...
        { none                      "None"                  false
        cs_facilitation            "CS fac"                true
        us_facilitation            "US fac"                true
        csus_facilitation         "CSUS fac"              true
        cs_fac_us_sup             "CS fac US sup"         true
        us_fac_cs_sup             "CS sup US fac"         true
        cs_suppression            "CS sup"                false
        us_suppression            "US sup"                false
        csus_suppression          "CSUS sup"              false
        };
    counts = [pieData{:,1}];
    labels = [pieData{:,2}];
    if withExplode
        explode = [pieData{:,3}];
    else
        explode = false(size(counts));
    end

elseif stim == "cs"
    [ none ...
        , cs_facilitation, us_facilitation ...
        , cs_suppression, us_suppression ...
        , csus_facilitation, csus_suppression ...
        , cs_fac_us_sup ...
        , us_fac_cs_sup ] = countChannelsForType(type, mouseCodes, stim);

    pieData = ...
        { none                      "None"                  false
        cs_facilitation            "CS fac"                true
        cs_suppression            "CS sup"                false
        };
    counts = [pieData{:,1}];
    labels = [pieData{:,2}];
    if withExplode
        explode = [pieData{:,3}];
    else
        explode = false(size(counts));
    end

else
    [ none ...
        , cs_facilitation, us_facilitation ...
        , cs_suppression, us_suppression ...
        , csus_facilitation, csus_suppression ...
        , cs_fac_us_sup ...
        , us_fac_cs_sup ] = countChannelsForType(type, mouseCodes, stim);

    pieData = ...
        { none                      "None"                  false
        us_facilitation            "US fac"                true
        us_suppression            "US sup"                false
        };
    counts = [pieData{:,1}];
    labels = [pieData{:,2}];
    if withExplode
        explode = [pieData{:,3}];
    else
        explode = false(size(counts));
    end
end

if contains(type, "MUT")
    if contains(type, "Shank2")
        colors = IkUtils.visualization.Params().red;
    elseif contains(type, "Tsc1KO")
        colors = IkUtils.visualization.Params().green;
    else
        colors = IkUtils.visualization.Params().blue;
    end
else
    if contains(type, "Shank2")
        colors = IkUtils.visualization.Params().darker_red;
    elseif contains(type, "Tsc1KO")
        colors = IkUtils.visualization.Params().darker_green;
    else
        colors = IkUtils.visualization.Params().darker_blue;
    end
end


% if contains(type, "MUT")
%     if contains(type, "Shank2")
%         colors = IkUtils.visualization.Params().red;
%     elseif contains(type, "Tsc1KO")
%         colors = IkUtils.visualization.Params().orange_red;
%     else
%         colors = IkUtils.visualization.Params().orange;
%     end
% else
%     if contains(type, "Shank2")
%         colors = IkUtils.visualization.Params().blue;
%     elseif contains(type, "Tsc1KO")
%         colors = IkUtils.visualization.Params().purple_blue;
%     else
%         colors = IkUtils.visualization.Params().purple;
%     end
% end

disp(pieData)
colormap(colors)

pie_chart = pie(counts, explode, labels);

% Remove labels for slices with zero values
text_handles = findobj(pie_chart, 'Type', 'text');
zero_values = counts == 0;
set(text_handles(zero_values), 'String', '');

% Adjust label positions to prevent overlap
for i = 1:numel(text_handles)
    text_handles = findobj(pie_chart, 'Type', 'text');
    % Check for overlap with other labels
    for j = i+1:numel(text_handles)
        if ~isempty(text_handles(i).String) && ~isempty(text_handles(j).String)
            while norm(text_handles(i).Position - text_handles(j).Position) < 0.26  % Adjust threshold as needed
                % Move the label to prevent overlap
                %                 text_handles(j).Position(1) = text_handles(j).Position(1) + 0.2;  % Adjust offset as needed
                if text_handles(i).Position(1) >= 0
                    if text_handles(i).Position(1) >= text_handles(j).Position(1)
                        text_handles(i).Position(1) = text_handles(i).Position(1) + 0.1;
                    else
                        text_handles(j).Position(1) = text_handles(j).Position(1) + 0.1;
                    end
                else
                    if text_handles(i).Position(1) >= text_handles(j).Position(1)
                        text_handles(i).Position(1) = text_handles(i).Position(1) - 0.1;  % Adjust offset as needed
                    else
                        text_handles(j).Position(1) = text_handles(j).Position(1) - 0.1;  % Adjust offset as needed
                    end
                end
            end
        end
    end
end
patch_handles = findobj(pie_chart, 'Type', 'patch');% Get the handles of patch objects (pie slices)

% add lines to labels
for i = 1:numel(text_handles)
    if ~isempty(text_handles(i).String)   % Check if the label is not empty
        label_position = text_handles(i).Position;% Get label position
        patch_position = get(patch_handles(i), 'Vertices');% Get corresponding patch (pie slice) position
        outer_midpoint = mean(patch_position([1, end], :));% Calculate the midpoint of the outer edge of the pie slice
        % Draw a leader line from the outer edge to the label
        line([outer_midpoint(1), label_position(1)], [outer_midpoint(2), label_position(2)], 'Color', 'k');
    end
end

uistack(patch_handles, 'top')

% Get current title position
current_position = get(get(gca, 'Title'), 'Position');

% Increase the vertical position
new_position = current_position + [0, 0.1, 0];

% Set the new title position
set(get(gca, 'Title'), 'Position', new_position);


nChannelsTotal = sum([none ...
    , cs_facilitation, us_facilitation ...
    , cs_suppression, us_suppression ...
    , csus_facilitation, csus_suppression ...
    , cs_fac_us_sup ...
    , us_fac_cs_sup ]);

title(sprintf("Modulation types of cells for %s", type))

try
    text ...
        ( fig.Children(1) ...
        , -1.7 ...
        , 1 ...
        , sprintf("N = %d", nChannelsTotal) ...
        , fontsize = 14 ...
        )
catch
end

end

function [ none ...
    , cs_facilitation, us_facilitation ...
    , cs_suppression, us_suppression ...
    , csus_facilitation, csus_suppression ...
    , cs_fac_us_sup ...
    , us_fac_cs_sup ] = countChannelsForType(type, mouseCodes, stim)

if nargin < 2
    mouseCodes = [defaultMice().code];
    mouseCodes = mouseCodes([defaultMice().type] == type);
end
if nargin < 3
    [noneVec ...
        , cs_facilitationVec, us_facilitationVec ...
        , cs_suppressionVec, us_suppressionVec ...
        , csus_facilitationVec, csus_suppressionVec ...
        , cs_fac_us_supVec ...
        , us_fac_cs_supVec] = cellfun ...
        ( @(mcode) countChannelsForMouse(mcode, type) ...
        , mouseCodes ...
        );

    none = sum(noneVec);
    cs_facilitation = sum(cs_facilitationVec);
    us_facilitation = sum(us_facilitationVec);
    cs_suppression = sum(cs_suppressionVec);
    us_suppression = sum(us_suppressionVec);
    csus_facilitation = sum(csus_facilitationVec);
    csus_suppression = sum(csus_suppressionVec);
    cs_fac_us_sup = sum(cs_fac_us_supVec);
    us_fac_cs_sup = sum(us_fac_cs_supVec);
elseif stim == "cs"

     [noneVec ...
        , cs_facilitationVec, us_facilitationVec ...
        , cs_suppressionVec, us_suppressionVec ...
        , csus_facilitationVec, csus_suppressionVec ...
        , cs_fac_us_supVec ...
        , us_fac_cs_supVec] = cellfun ...
        ( @(mcode) countChannelsForMouse(mcode, type) ...
        , mouseCodes ...
        );

    none = sum(noneVec) + sum(us_facilitationVec) + sum(us_suppressionVec);
    cs_facilitation = sum(cs_facilitationVec) + sum(csus_facilitationVec) + sum(cs_fac_us_supVec);
    cs_suppression = sum(cs_suppressionVec) + sum(csus_suppressionVec) + sum(us_fac_cs_supVec);    
    us_facilitation =  [];
    us_suppression =  [];
    csus_facilitation = [];
    csus_suppression = [];
    cs_fac_us_sup = [];
    us_fac_cs_sup = [];

    
else
  [noneVec ...
        , cs_facilitationVec, us_facilitationVec ...
        , cs_suppressionVec, us_suppressionVec ...
        , csus_facilitationVec, csus_suppressionVec ...
        , cs_fac_us_supVec ...
        , us_fac_cs_supVec] = cellfun ...
        ( @(mcode) countChannelsForMouse(mcode, type) ...
        , mouseCodes ...
        );

    none = sum(noneVec) + sum(cs_facilitationVec) + sum(cs_suppressionVec);
    cs_facilitation = [];
    cs_suppression = [];    
    us_facilitation = sum(us_facilitationVec) + sum(csus_facilitationVec) + sum(us_fac_cs_supVec);
    us_suppression = sum(us_suppressionVec) + sum(csus_suppressionVec) + sum(cs_fac_us_supVec);
    csus_facilitation = [];
    csus_suppression = [];
    cs_fac_us_sup = [];
    us_fac_cs_sup = [];

end




end

function [ none ...
    , cs_facilitation, us_facilitation ...
    , cs_suppression, us_suppression ...
    , csus_facilitation, csus_suppression ...
    , cs_fac_us_sup ...
    , us_fac_cs_sup] = countChannelsForMouse(mcode, type)

% sessionIdcs = splitSessionsByType(mcode).(type);

% allModulationMasks = loadSimpleChannelModulationMasks(mcode);
% allSimpleChannelMasks = loadSimpleMasks(mcode);
%
% modulationMasks = allModulationMasks(sessionIdcs);
% SimpleMasks = allSimpleChannelMasks(sessionIdcs);

modulationMasks = loadSimpleChannelModulationMasks(mcode);
simpleMasks = loadSimpleMasks(mcode);
if isempty(modulationMasks)
    modulationMasks = {[]};
end
if isempty(simpleMasks)
    simpleMasks = {[]};
end

[noneVec ...
    , cs_facilitationVec, us_facilitationVec ...
    , cs_suppressionVec, us_suppressionVec ...
    , csus_facilitationVec, csus_suppressionVec ...
    , cs_fac_us_supVec ...
    , us_fac_cs_supVec] = cellfun ...
    ( @(modMask, compMask) countChannelsForSession(modMask(compMask)) ...
    , modulationMasks ...
    , simpleMasks ...
    );

none = sum(noneVec);
cs_facilitation = sum(cs_facilitationVec);
us_facilitation = sum(us_facilitationVec);
cs_suppression = sum(cs_suppressionVec);
us_suppression = sum(us_suppressionVec);
csus_facilitation = sum(csus_facilitationVec);
csus_suppression = sum(csus_suppressionVec);
cs_fac_us_sup = sum(cs_fac_us_supVec);
us_fac_cs_sup = sum(us_fac_cs_supVec);
end

function [ none ...
    , cs_facilitation, us_facilitation ...
    , cs_suppression, us_suppression ...
    , csus_facilitation, csus_suppression ...
    , cs_fac_us_sup ...
    , us_fac_cs_sup] = countChannelsForSession(modulationMasks)

%          keyboard

if isempty(modulationMasks)
    none = 0;
    cs_facilitation = 0;
    us_facilitation = 0;
    cs_suppression = 0;
    us_suppression = 0;
    csus_facilitation = 0;
    csus_suppression = 0;
    cs_fac_us_sup = 0;
    us_fac_cs_sup = 0;
    return
end

[noneVec ...
    , cs_facilitationVec, us_facilitationVec ...
    , cs_suppressionVec, us_suppressionVec ...
    , csus_facilitationVec, csus_suppressionVec ...
    , cs_fac_us_supVec ...
    , us_fac_cs_supVec] = cellfun ...
    ( @parseChannelMask ...
    , modulationMasks ...
    );

none = sum(noneVec);
cs_facilitation = sum(cs_facilitationVec);
us_facilitation = sum(us_facilitationVec);
cs_suppression = sum(cs_suppressionVec);
us_suppression = sum(us_suppressionVec);
csus_facilitation = sum(csus_facilitationVec);
csus_suppression = sum(csus_suppressionVec);
cs_fac_us_sup = sum(cs_fac_us_supVec);
us_fac_cs_sup = sum(us_fac_cs_supVec);

end

function [ none ...
    , cs_facilitation, us_facilitation ...
    , cs_suppression, us_suppression ...
    , csus_facilitation, csus_suppression ...
    , cs_fac_us_sup ...
    , us_fac_cs_sup] = parseChannelMask(mask)

none = 0;
cs_facilitation = 0;
us_facilitation = 0;
cs_suppression = 0;
us_suppression = 0;
csus_facilitation = 0;
csus_suppression = 0;
cs_fac_us_sup = 0;
us_fac_cs_sup = 0;

if isempty(mask)
    return
elseif mask.cs_facilitation
    if mask.us_facilitation
        csus_facilitation = 1;
    elseif mask.us_suppression
        cs_fac_us_sup = 1;
    else
        cs_facilitation = 1;
    end
elseif mask.us_facilitation
    if mask.cs_suppression
        us_fac_cs_sup = 1;
    else
        us_facilitation = 1;
    end
elseif mask.cs_suppression
    if mask.us_suppression
        csus_suppression = 1;
    else
        cs_suppression = 1;
    end
elseif mask.us_suppression
    us_suppression = 1;
else
    none = 1;
end

end

