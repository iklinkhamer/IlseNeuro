%% IK 29-5-24
function masks = simpleCurationMouseWrapperNoUser()

mouseCodes = arrayfun ...
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

for m = 1:length(mouseCodes)
    mcode = mouseCodes(m);
    mname = mouseNames(m);
    spikeSortedEphysData = getData(mcode);

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
            fprintf("There are no neurons for %s, %s.\n", mname, mcode)
            masks = [];
            continue
        end
    catch
        fprintf("There are no neurons for %s, %s.\n", mname, mcode)
        masks = [];
        continue
    end


    candidateMasks = computeChannelCandidates(spikeSortedEphysData);

    nSessions = length(spikeSortedEphysData);
    allSessionIdcs = 1:nSessions;

    for s = 1:nSessions
        for n = 1:length(spikeSortedEphysData(s).neuron)
            spikeSortedEphysData(s).neuron(n).RasterXY_us = spikeSortedEphysData(s).neuron(n).RasterXY_cs;
        end
    end

    axs = struct ...
        ( neuronAxs = IkUtils.initPlots([1 3]) ...
        );


    %     prevCuration = loadSimpleCurationResultRaw(mcode, onlyLatest = true);
    prevCuration = [];

    if numel(prevCuration) ~= nSessions
        if numel(prevCuration) > 1
            warning ...
                ( "No. sessions in previous curation result does not match the " ...
                + "no. sessions found for %s:\n" ...
                + "\tprevious: %d\n" ...
                + "\tcurrent: %d\n" ...
                , mcode ...
                , numel(prevCuration) ...
                , nSessions ...
                )
            keyboard
        else
            prevCuration = @(~) [];
        end
    end

    masks = arrayfun ...
        ( @(sessionIdx) ...
        simpleCurationSessionWrapperNoUser ...
        ( axs ...
        , mname ...
        , sessionIdx ...
        , spikeSortedEphysData(sessionIdx)...
        , candidateMasks(sessionIdx).simpleMask ...
        , prevCuration(sessionIdx) ...
        ) ...
        , allSessionIdcs ...
        );

    try
        saveSimpleCurationResultNoUser(mcode, masks)
    catch
        keyboard
    end
end
end
