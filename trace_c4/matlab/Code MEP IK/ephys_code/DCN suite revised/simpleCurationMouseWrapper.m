%% I.K. 1-6-24
% Curate complex channels of all sessions for the selected mouse
function masks = simpleCurationMouseWrapper(mname)
        
    if nargin < 1
        
        fprintf("\n-------------------------------------")
        fprintf("\nWhich mouse would you like to curate?")
        fprintf("\n-------------------------------------\n")
        
        mouseCodes = arrayfun ...          % IK change
            ( @(mouse) string(mouse.code) ...
            , defaultMice() ...
            );
        mouseNames = arrayfun ...
            ( @(mouse) string(mouse.name) ...
            , defaultMice() ...
            );

        ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
        ephysMiceMask(end-1:end) = true;
        mouseNames = mouseNames(ephysMiceMask);
        mouseCodes = mouseCodes(ephysMiceMask);
        %mouseNames = [mouseNames, "All of the above"];
        mouseNames = mouseCodes + "     " + mouseNames;
        
        [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
        mcode = mouseCodes(mname_idx);
        %mname = mouseNames(mname_idx);
    end
    
%     spikeSortedEphysData = loadAnalyzedEphysDataForMouse(mname);
    spikeSortedEphysData = getData(mcode); % IK
    
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

    
    candidateMasks = computeChannelCandidates(spikeSortedEphysData);
    
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
