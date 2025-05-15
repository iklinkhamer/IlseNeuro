%% I.K. 1-6-24
% Curate complex channels of all sessions for the selected mouse
function masks = simpleCurationMouseWrapper_CSG(mname, kwargs)
arguments
    mname = '';
    kwargs.directory = fullfile(Env.getBayesLabUserRoot, "/ContextMouseExperiments/Ilse/CSG/");
    kwargs.data_directory = fullfile(Env.getBayesLabUserRoot, "/ContextMouseExperiments/Analysis/BEPAK/analysis_directories/")
    kwargs.batchSelection = 1;
end

if isempty(mname)

    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to curate?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...          % IK change
        ( @(mouse) string(mouse.code) ...
        , defaultMice_CSG() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice_CSG() ...
        );

    ephysMiceMask = ~cellfun(@isempty,{defaultMice_CSG().ephysdates});
    ephysMiceMask(end-1:end) = true;
    mouseNames = mouseNames(ephysMiceMask);
    mouseCodes = mouseCodes(ephysMiceMask);
    %mouseNames = [mouseNames, "All of the above"];
    mouseNames = mouseCodes + "     " + mouseNames;

    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    mcode = mouseCodes(mname_idx);
        mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice_CSG() ...
        );
    mname = mouseNames(mname_idx);
end



getSessionsLater = 1;
%i= 0;

ephys_session_dir = dir(fullfile(kwargs.directory, sprintf('StructEphysData_%s*.mat', mname)));
ephys_session_files = string({ephys_session_dir.name});
% ephysDataVec_ = [];
% for f = 1:length(ephys_session_files)
%     csg_data = load(fullfile(kwargs.dataFolder, ephys_session_files(f)));
%     neuron_struct = csg_data.unit.neuron;
%     ephysDataVec_= [ephysDataVec_, neuron_struct];
% end

sessions = 1:length(ephys_session_files);

mdir = dir(fullfile(kwargs.data_directory, mname));

    folder_names = {mdir.name};
    numeric_suffix = cellfun(@(n) regexp(n, '\d+$', 'match', 'once'), folder_names, 'UniformOutput', false);
    numeric_vals = cellfun(@(s) str2double(s), numeric_suffix);
    if getSessionsLater
        sessions = numeric_vals(~isnan(numeric_vals));
    end
    valid_indices = [mdir.isdir] & ~ismember(folder_names, {'.', '..'}) & ismember(numeric_vals, sessions);
    mfolders = string({mdir(valid_indices).name});

for s=1:length(sessions)
    folder = mfolders(s);
    % session = sessions(:)'
   % session = sessions(s)';
    %i = i+1;
    %if i ~=2
    %     continue
    %end
    % [fig, ax] = viewRaster(neuron, "cs", trialSubset = 1:80, ax = rasterAxes(1)
    %   


    try % AnalyzedEphys might not exist
        csg_data = load(fullfile(kwargs.directory, ephys_session_files(s)));
        units = csg_data.unit;
        units.neuron = units.neuron(kwargs.batchSelection:end);
        % neurons = csg_data.unit.neuron;
        % stimtimes = csg_data.unit.stimtimes;
    catch       
        sprintf("Loading units failed for session %d", s)
        continue
    end
%     spikeSortedEphysData = loadAnalyzedEphysDataForMouse(mname);
    % spikeSortedEphysData = getData(mcode); % IK
    spikeSortedEphysData(s) = units;

end
    
    try 
        if isempty(spikeSortedEphysData(2).neuron)
            empty_s2 = 1;
        else
            empty_s2 = 0;
        end
    catch
        empty_s2 = 1;
    end
try
    if isempty(spikeSortedEphysData(1).neuron) && empty_s2
        fprintf("There are no neurons for this mouse.")
        masks = [];
        return
    end
catch
    fprintf("There are no neurons for this mouse.")
        masks = [];
        return
end

    
    candidateMasks = computeChannelCandidates_CSG(spikeSortedEphysData);
    
    %sessionIdcs = splitSessionsByType(mname);
    %allSessionIdcs = [sessionIdcs.delta sessionIdcs.uniform];
    %nSessions = numel(allSessionIdcs);%endSession - startSession + 1;
    nSessions = length(spikeSortedEphysData); % IK
    allSessionIdcs = 1:nSessions;

       for s = 1:nSessions % IK change
           for n = 1:length(spikeSortedEphysData(s).neuron)
                spikeSortedEphysData(s).neuron(n).RasterXY_us = spikeSortedEphysData(s).neuron(n).RasterXY_cs;
           end
       end
    
    axs = struct ...
        ( neuronAxs = IkUtils.initPlots([1 3]) ... % IK change
        );
    

    prevCuration = loadSimpleCurationResultRaw(mcode, onlyLatest = true);
    
    if numel(prevCuration) ~= nSessions
        if numel(prevCuration) > 1 % IK change
            warning ...
                ( "No. sessions in previous curation result does not match the " ...
                + "no. sessions found for %s:\n" ...
                + "\tprevious: %d\n" ...
                + "\tcurrent: %d\n" ...
                , mname ...
                , numel(prevCuration) ...
                , nSessions ...
                )
            keyboard
        else
            fprintf("\nNo previous curation results found for %s.\n", mname)
            prevCuration = @(~) [];
        end
    end
        
    masks = arrayfun ...
        ( @(sessionIdx) ...
            simpleCurationSessionWrapper ...
            ( axs ...
            , mname ...
            , sessionIdx ...
            , spikeSortedEphysData(sessionIdx)...
            , candidateMasks(sessionIdx).simpleMask ...
            , prevCuration(sessionIdx) ...
            ) ...
        , allSessionIdcs ...     
        );

    saveSimpleCurationResult(mcode, masks) 
    
  
    
end
