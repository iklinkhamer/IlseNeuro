%% Ilse Klinkhamer 28-02-2025
function inspectTrials(mouseName, kwargs)
    arguments
        mouseName string = "Greene";
        kwargs.directory string = fullfile(Env.getBayesLabUserRoot, "TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder/", mouseName);
    end

    % Adjust directory if mouseName contains "Reserve"
    if contains(mouseName, "Reserve")
        kwargs.directory = strrep(kwargs.directory, "MainFolder", "ReserveFolder");
    end
    single_pattern = [1,1,1,1,1];

    mouse = Subject(mouseName);
    sessions = mouse.collectSessions(sessionkinds="EPHYS");
    
    fprintf("Found %d sessions for %s\n\n", numel(sessions), mouseName)

    for session = sessions(:)'
        timestamp_sess = session.timestampIdStr;
        session_name = sprintf("%s_%s", mouseName, timestamp_sess);
        filename = sprintf("EphysSession_%s.mat", session_name);
        
        filePath = fullfile(kwargs.directory, session_name, filename);

        if isfile(filePath)
            data = load(filePath);
            trials = data.param.exp.pairindex
            if any(trials ~= single_pattern)
                firstWideSession = sprintf("%s_%s", mouseName, session.timestampIdStr);
                outputFile = fullfile(kwargs.directory, "firstWideSession.txt");

                % Save session info to a text file
                fid = fopen(outputFile, 'w');
                if fid ~= -1
                    fprintf(fid, "%s\n", firstWideSession);
                    fclose(fid);
                else
                    warning("Could not open file for writing: %s", outputFile);
                end
                return
            end
        else
            warning("File not found: %s", filePath);
        end
    end
end