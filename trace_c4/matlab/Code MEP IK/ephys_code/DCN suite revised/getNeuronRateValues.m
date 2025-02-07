function neuron_rates = getNeuronRateValues(sessionsData)

nNeuronsPerSession = getNeuronsPerSession(sessionsData);

%neuron_base_rates = zeros(nNeurons,1);
%neuron_max_rates = zeros(nNeurons,1);

neuron_rates = struct;

nSessions = length(sessionsData);


for s = 1:nSessions
    
    % TODO: Julius: splice out the `index` from the `getConditions` function
    [nConditions, index, ~, ~] = getConditions(sessionsData(s));
    
    for n = 1 : nNeuronsPerSession(s)
        
        %n_base_rate = mean(mean(data(s).neuron(n).psth_cs_reset( : , index.base_idx_cs)));
        %neuron_base_rates((nNeurons_prev_sessions(s) + n),1) = n_base_rate;
        neuron_rates(s).base_rates(n) = ...
            mean ...
                ( mean ...
                    ( sessionsData(s).neuron(n).psth_cs_reset(:, index.base_idx_cs) ...
                    ) ...
                );
        
        %n_max_rate = median(max(data(s).neuron(n).psth_cs_reset(2 : nConditions,index.maxmaxrange_cs),[],2)) ;
        %neuron_max_rates((nNeurons_prev_sessions(s) + n) ,1) = n_max_rate;
        neuron_rates(s).max_rates(n) = ...
            median ...
                ( max ...
                    ( sessionsData(s).neuron(n).psth_cs_reset(2:nConditions, index.maxmaxrange_cs) ...
                    , [] ...
                    , 2 ...
                    ) ...
                );
    end
    
    
end

end