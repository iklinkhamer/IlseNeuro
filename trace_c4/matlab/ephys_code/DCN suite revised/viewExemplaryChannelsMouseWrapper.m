%% IK 28-5-24
function viewExemplaryChannelsMouseWrapper()
app = getUImain(scriptname=mfilename);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    appData = evalin('base', 'appData');

    includeUSaligned = appData.includeUSaligned;

    if appData.enumeration == "Single mouse"
        mcode = appData.mouseCode;
    else
        disp("Pick a different enumeration")
        continue
    end

    simpleCuration = loadSimpleCurationResultSanitized(mcode); % IK removed complexCuration

    ephysData = loadAnalyzedEphysDataForMouse(mcode);
    nSessions = numel(simpleCuration); % IK change

    if ~includeUSaligned
        axs = IkUtils.initPlots([1 3]); % IK change
    else
        axs = IkUtils.initPlots([2 3]); % IK change
    end

    arrayfun ...
        ( @(s) viewDataSessionWrapper ...
        ( axs ...
        , ephysData(s) ...
        , simpleCuration(s).exemplary ... % IK removed complexCuration
        , s ...
        , mcode ...
        ) ...
        , 1:nSessions ...
        )
    disp("End of loop, waiting for UI input.")
end

end

% IK removed myOr function