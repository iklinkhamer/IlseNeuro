%% IK 29-5-24
function maskStruct = simpleCurationSessionWrapperNoUser ...
    ( axs ...
    , mname ...
    , sessionIdx ...
    , sessionData ...
    , simpleCandidateMask ...
    , previousCuration ...
    )

nChannels = numel(simpleCandidateMask);

p = IkUtils.getParams();

if ~any(simpleCandidateMask)
    warning ...
        ( "No simple spike candidates found in session %d for %s" ...
        , sessionIdx ...
        , mname ...
        );
    
    if isempty(sessionData.neuron) 
        maskStruct = struct ...
            ( simple = simpleCandidateMask ...
            , exemplary = false(1, nChannels) ... 
            , modulation = {cell(1,nChannels)} ...
            );

        for i = 1 : length(maskStruct.modulation)
            if maskStruct.simple(i) == 0
                maskStruct.modulation{i} = struct ...
                    ( cs_facilitation = false ...
                    , us_facilitation = false ...
                    , cs_suppression = false ...
                    , us_suppression = false ...
                    );
            end
        end
        return
    end
end

if isempty(previousCuration)
    masks = struct ...
        ( simple = simpleCandidateMask ...
        , exemplary = false(1, nChannels) ... 
        , modulation = {cell(1,nChannels)} ... 
        );
    for i = 1 : length(masks.modulation)
        if masks.simple(i) == 0
            masks.modulation{i} = struct ...
                ( cs_facilitation = false ...
                , us_facilitation = false ...
                , cs_suppression = false ...
                , us_suppression = false ...
                );
        end
    end
else
    masks = previousCuration;
end


simpleCandidateMask = masks.simple;

simpleCandidates = find(simpleCandidateMask);

if isempty(simpleCandidates)
    uiState = UiState ...
        ( channelIdx = 0 ...
        , masks = masks ...
        , reachedBottom = true ...
        , reachedTop = false ...
        );
    maskStruct = uiState.masks;
    return
end

choiceMap = ...
    { "toggle simple selection"     , @toggleSimple ...
    ; "override modulation"         , @promptModulationMask ...
    ; "mark as exemplary"           , @toggleExemplaryNeuron ...
    ; "debug"                       , @debug ...
    ; "prev"                        , mkPrevFn(simpleCandidateMask) ...
    ; "next"                        , mkNextFn(simpleCandidateMask) ...
    };
quit = false;

uiState = UiState ...
    ( channelIdx = simpleCandidates(1) ... 
    , masks = masks ...
    , reachedBottom = false ...
    , reachedTop = false ...
    );

while ~quit

    choice = "";

    while not( uiState.reachedBottom || uiState.reachedTop )% || uiState.exit )

        %% Compute histogram statistics and estimate modulation
        modulationRanges = IkUtils.getParams().sspkRanges;

        histStats = struct ...
            ( cs = computePsthStats ...
            ( sessionData.neuron(uiState.channelIdx).RasterXY_cs_filtered ...
            , modulationRanges.cs ...
            , p.psthRanges.cs_full ...
            ) ... 
            , us = computePsthStats ...
            ( sessionData.neuron(uiState.channelIdx).RasterXY_cs_filtered ... 
            , modulationRanges.us ...
            , p.psthRanges.cs_full ... 
            ) ...
            );

        if isempty(uiState.masks.modulation{uiState.channelIdx})
            uiState.masks.modulation{uiState.channelIdx} = struct ...
                ( cs_facilitation = histStats.cs.correctedRate > histStats.cs.baserate + 8 ... 
                , cs_suppression = histStats.cs.correctedRate < histStats.cs.baserate - 8 ...
                , us_facilitation = histStats.us.correctedRate > histStats.us.baserate + 8 ... 
                , us_suppression = histStats.us.correctedRate < histStats.us.baserate - 8 ... 
                );
        end        

        %% Check for and visualize previous curation

        if ~uiState.masks.simple(uiState.channelIdx)
            arrayfun ...
                ( @(ax) fadeAxes(ax, color = [0.7048 0.5418 0.5418]) ...
                , axs.neuronAxs ...
                )
            fprintf("\n\nThis channel is currently not in the simple channel selection.\n\n")
        end

        if uiState.masks.exemplary(uiState.channelIdx)
            axs.neuronAxs(1).Parent.Color = getColors().exemplary;
        else
            axs.neuronAxs(1).Parent.Color = [.98 .98 .98];
        end

        %% CURATE automatically

        fn = choiceMap{6, 2};

        try
            uiState = fn(uiState);
        catch err
            warning ...
                ( "Caught error during execution of `%s`:\n%s\n" ...
                + "\terror: %s\n" ...
                + "\tfunction: %s\n" ...
                + "\tline: %d\n" ...
                , choice ...
                , fprintf(string(err.message)) ... % IK change
                , sprintf("%s: %s", err.identifier, err.message) ...
                , err.stack(1).name ...
                , err.stack(1).line ...
                )

            fprintf("\nStarting debug mode.\n\n")
            keyboard
        end

    end

    if uiState.reachedBottom
        warning("Already at first neuron of session")
        uiState.reachedBottom = false;
    end
    if uiState.reachedTop
        quit = true;
    end

end

maskStruct = uiState.masks;

end

function uiState = toggleExemplaryNeuron(uiState)
arguments
    uiState(1,1) UiState
end

idx = uiState.channelIdx;
exemplaryStatus = uiState.masks.exemplary(idx);
uiState.masks.exemplary(idx) = not(exemplaryStatus);

if uiState.masks.exemplary(idx)
    uiState.masks.simple(idx) = true;
end

end

function uiState = promptModulationMask(uiState)

modulationCombinations = ...
    { "none"            , [] ...
    ; "CS facilitation only"         , ["cs_facilitation"] ...
    ; "US facilitation only"         , ["us_facilitation"] ...
    ; "CS suppression only"          , ["cs_suppression"] ...
    ; "US suppression only"          , ["us_suppression"] ...
    ; "CS and US facilitation"       , ["cs_facilitation", "us_facilitation"] ... % IK change, removed prior
    ; "CS facilitation and US suppression"  , ["cs_facilitation", "us_suppression"] ...
    ; "CS suppression and US facilitation"  , ["cs_suppression", "us_facilitation"] ...
    ; "CS and US suppression"               , ["cs_suppression", "us_suppression"] ...
    };

[~, idx] = IkUtils.do_prompt_select_option([modulationCombinations{:,1}]);

newMask = modulationStruct(modulationCombinations{idx, 2});

uiState.masks.modulation{uiState.channelIdx} = newMask;

    function result = modulationStruct(events)

        result = struct ...
            ( cs_facilitation = false ... % IK change, removed prior
            , us_facilitation = false ...
            , cs_suppression = false ...
            , us_suppression = false ...
            );

        for event = events(:)'
            result.(event) = true;
        end

    end

end

function uiState = toggleSimple(uiState)
arguments
    uiState(1,1) UiState
end

idx = uiState.channelIdx;
selectionStatus = uiState.masks.simple(idx);
uiState.masks.simple(idx) = not(selectionStatus);

if not(uiState.masks.simple(idx))
    uiState.masks.exemplary(idx) = false;
end
end

function uiState = debug(uiState)
keyboard
end

%% Navigation Utilities

function bool = hasPrev(candidateMask, idx)
bool = idx > find(candidateMask, 1);
end

function bool = hasNext(candidateMask, idx)
bool = idx < find(candidateMask, 1, 'last');
end

function fn = mkPrevFn(candidateMask)

    function uiState = prevChannel(uiState)
        currentIdx = uiState.channelIdx;
        if hasPrev(candidateMask, currentIdx)
            prevIdx = find(candidateMask(1:(currentIdx-1)), 1, 'last');
        else
            prevIdx = currentIdx;
            uiState.reachedBottom = true;
        end
        uiState.channelIdx = prevIdx;
    end

fn = @prevChannel;
end

function fn = mkNextFn(candidateMask)

    function uiState = nextChannel(uiState)
        currentIdx = uiState.channelIdx;
        if hasNext(candidateMask, currentIdx)
            nextIdx = currentIdx + find(candidateMask(currentIdx+1:end, 1), 1);
        else
            nextIdx = currentIdx;
            uiState.reachedTop = true;
        end
        uiState.channelIdx = nextIdx;
    end

fn = @nextChannel;
end