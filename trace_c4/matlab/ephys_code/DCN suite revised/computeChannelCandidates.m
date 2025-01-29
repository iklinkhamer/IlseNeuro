%% I.K. 1-6-24
function candidateMask = computeChannelCandidates(data)

    clusterMask = computeNeuronClusters(data);

    contamination = contaminationFilter(data); % IK note

    numValidTrialsAfterThres = arrayfun(@(session) computeValidTrials(session), data, UniformOutput=false);

    candidateMask = struct;

    for s = 1 : length(data)
        candidateMask(s).simpleMask = clusterMask(s).simpleMask(:);

        try
        if isempty([contamination(s).neuron.percentage])
            continue
        end
        catch
            continue
        end

        contaminationMask = [contamination(s).neuron.percentage] <= 10;
        trialMask = numValidTrialsAfterThres{s} > IkUtils.getParams().validTrialThres;
        candidateMask(s).simpleMask = clusterMask(s).simpleMask(:) ...
                                     & contaminationMask(:) & trialMask(:);                             
    end

end
function contamination = contaminationFilter(sessionsData)
    contamination = arrayfun ...
        ( @(session) computeContaminationForSession(session.neuron) ...
        , sessionsData ...
        );    
end

function contamination = computeContaminationForSession(neurons)
    
    contamination.neuron = arrayfun ...
        ( @computeContamination ...
        , neurons ...
        );
    
end
function contamination = computeContamination(channelData)
    
% %     spikeTimes = channelData.st; % IK change
     spikeTimes = channelData.spiketimes; % IK change

    contamination.spike_intervals = diff(spikeTimes) * 1000;

    shortIntervals = sum(contamination.spike_intervals < 2); % IK change from 5 to 1. 18/12/23
    
    totalIntervals = numel(contamination.spike_intervals);
    
    contamination.percentage = shortIntervals / totalIntervals * 100;
    
    fprintf("Neuron: %d Contamination percentage: %2f%%\n", channelData.irc_trueID, contamination.percentage);

end


function validTrials = computeValidTrials(sessionData)
    validTrials = [];
    for n = 1 : length(sessionData.neuron)
        validTrials(n) = numel(unique(sessionData.neuron(n).RasterXY_cs_filtered(2,1:3:end)));
        fprintf("Neuron: %d Number of valid trials: %d\n", n, validTrials(n))
    end
end