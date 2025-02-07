%% I.K. 1-6-24
function mask = computeNeuronClusters(sessionsData)

neuron_rates = getNeuronRateValues(sessionsData);

nCumulativeNeuronsPrevSessions = [0 getCumulativeNeuronsPerSession(sessionsData)];

nSessions = length(sessionsData);

neuron_base_rates = [];

for s = 1 : nSessions
    neuron_base_rates = [neuron_base_rates, neuron_rates(s).base_rates];
end

for i = 1 : length(neuron_base_rates)
    fprintf("neuron: %d base rate: %2f\n", i, neuron_base_rates(i))
end

SS_base_filter_mask = neuron_base_rates >= IkUtils.getParams().Sspk_base_filter;

simple_mask_all_sessions = SS_base_filter_mask(:);
mask = struct;

for s = 1 : nSessions
    
    firstNeuron = nCumulativeNeuronsPrevSessions(s)+1;
    lastNeuron = nCumulativeNeuronsPrevSessions(s+1);
    
    mask(s).simpleMask = simple_mask_all_sessions(firstNeuron:lastNeuron);
    
end
 
end


function neuron_rates = getNeuronRateValues(sessionsData)

nNeuronsPerSession = getNeuronsPerSession(sessionsData);

neuron_rates = struct;

nSessions = length(sessionsData);


for s = 1:nSessions
    
    if ~isempty(sessionsData(s).neuron)
        [nConditions, index, ~, ~] = getConditions(sessionsData(s));
    end
    
    for n = 1 : nNeuronsPerSession(s)
        neuron_rates(s).base_rates(n) = ...
            nanmean ...
                ( nanmean ...
                    ( sessionsData(s).neuron(n).psth_cs_reset_filtered(:, index.base_idx_cs) ... % IK change 3-4-24
                    ) ...
                );
        neuron_rates(s).max_rates(n) = ...
            nanmedian ...
                ( nanmax ...
                    ( sessionsData(s).neuron(n).psth_cs_reset_filtered(2:nConditions, index.maxmaxrange_cs) ...
                    , [] ...
                    , 2 ...
                    ) ...
                );
    end
    
    
end

end

function nCumPerSess = getCumulativeNeuronsPerSession(sessionsData)
    
    nCumPerSess = cumsum(getNeuronsPerSession(sessionsData));
    
end

function nNeuronsPerSession = getNeuronsPerSession(sessionsData)

for i = 1:length(sessionsData) 
   nNeuronsPerSession(i) = length(sessionsData(i).neuron);
end

end