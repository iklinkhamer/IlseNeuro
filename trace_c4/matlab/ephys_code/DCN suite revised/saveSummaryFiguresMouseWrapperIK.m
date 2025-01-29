%% IK 28-5-24
function saveSummaryFiguresMouseWrapperIK()
close all
app = getUImain(scriptname=mfilename);
disp("Waiting for UI input...")

loadPrevData = true;

while true    
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    disp("Calculating...")
    appData = evalin('base', 'appData');
    gene = appData.gene;
    loc = appData.loc;
    
    mouseCodes = appData.mousecodes;
    mouseTypes = appData.mousetypes;

    mouseCodesDefaultMice = arrayfun ...          % IK change
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
    ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
    ephysMiceMask(end-1:end) = true;
    mouseCodeEphysIdcs = contains(mouseCodes, mouseCodesDefaultMice(ephysMiceMask));

    mouseCodes = mouseCodes(mouseCodeEphysIdcs);
    mouseTypes = mouseTypes(mouseCodeEphysIdcs);
    WTmice = mouseCodes(contains(mouseTypes, gene+loc+"WT"));
    MUTmice = mouseCodes(contains(mouseTypes, gene+loc+"MUT"));

    mice = {WTmice, MUTmice};

    P = IkUtils.getParams();
    outputRoot = fullfile ...
        (fullfile(P.figPath, 'boxPlotsNew'));

    if ~exist(outputRoot, 'dir')
        mkdir(outputRoot)
    end

    plotSummaryFiguresMUTvsWT_IK ...
        ( mice ...
        , gene ...
        , loc ...
        , outputRoot...
        , loadPrevData...
        );

    disp("Loop done. Waiting for new UI input...")
end

end