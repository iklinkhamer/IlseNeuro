%% IK 15-3-2024
function overviewTracesSingleMouse(mInput)
savefigures = true;
startMouse = "MI24.00245.04";
day = 10;

if nargin < 1
    mouseCodes = IkUtils.promptMouseNames();
else
    mouseCodes = IkUtils.promptMouseNames(mInput);
end
non_valid_idcs = find(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT");
mouseCodes(non_valid_idcs) = [];

if isstring(startMouse) || ischar(startMouse)
    startMouse = find(mouseCodes == startMouse);
else
    if ~isempty(find(startMouse>non_valid_idcs, 1, 'last' ))
        startMouse = startMouse - find(startMouse>non_valid_idcs, 1, 'last' );
    end
end
if length(mouseCodes) == 1
    startMouse = 1;
end
for m = startMouse:length(mouseCodes)
    close all    
    dataStruct = calcMouseForOverviewTraces(mouseCodes(m)); % calculations    
    plotEyeTracesFigures(dataStruct, mouseCodes(m), savefigures, day) % Final figures
    plotCRoverview(dataStruct, mouseCodes(m), savefigures, day)
end
end

function plotEyeTracesFigures(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotEyeTracesFigures')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end

eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotEyelidTraces(Mouse, mcode, eyeblinkTracesOverview, [2,1,1], Day)
plotEyeTracesOverTime(Mouse, mcode, eyeblinkTracesOverview, [2,1,2], Day)

if savefigures
    fname = sprintf('eyeblinkTracesOverview_%s.eps', mcode);
    file = fullfile(P.figPath, "Training overview eps",fname);
    saveas(eyeblinkTracesOverview,file);
    fname2 = sprintf('eyeblinkTracesOverview_%s.png', mcode);
    file = fullfile(P.figPath, "Training overview",fname2);
    saveas(eyeblinkTracesOverview,file);
end
end

function plotCRoverview(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotCRoverview')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure(Mouse, mcode, CROverview, [2,2,1], day = Day);
plotCRperc(Mouse, mcode, CROverview, [2,2,2], day = Day);
plotOffset(Mouse, mcode, CROverview, [2,2,3], Day);

if savefigures
    fname = sprintf('CRoverview_%s.eps', mcode);
    file = fullfile(P.figPath, "Training overview eps",fname);
    saveas(CROverview,file);
    fname2 = sprintf('CRoverview_%s.png', mcode);
    file = fullfile(P.figPath, "Training overview",fname2);
    saveas(CROverview,file);
end
end