function nCumPerSess = getCumulativeNeuronsPerSession(sessionsData)
    
    nCumPerSess = cumsum(getNeuronsPerSession(sessionsData));
    
end

function nNeuronsPerSession = getNeuronsPerSession(sessionsData)

    nNeuronsPerSession = arrayfun ...
        ( @(session) size(session.AveFiring_cs, 1) ...
        , sessionsData ...
        );

end