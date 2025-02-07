
function switchTypes = detectSwitchSessions(sessionsData)

   

    function distType = getDistType(session)
        distType_ = sessionType(session);
        distType = string(distType_);
    end

    distTypes = arrayfun(@getDistType, sessionsData);

    nonSwitchVal = struct(mask = false, type = "");

    function switchType = getSwitchType(prevDist, currentDist)
        if strcmp(prevDist, currentDist)
            switchType = nonSwitchVal;
            return
        end

        prevDistChar = char(prevDist);
        curDistChar = char(currentDist);

        switchType = struct ...
                ( mask = true ...
                , type = strcat(prevDistChar(1), curDistChar(1)) ...
                );
    end

    switchTypes_ = arrayfun ...
        ( @getSwitchType ...
        , distTypes(1:end-1) ... % 'previous' sessions
        , distTypes(2:end) ... % 'current' sessions
        );

    switchTypes = ...
        [ nonSwitchVal ... % first session can't be a switch by definition
        , switchTypes_ ...
        ];

end