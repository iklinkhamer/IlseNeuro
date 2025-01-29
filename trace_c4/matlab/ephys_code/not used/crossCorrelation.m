function n_simple = crossCorrelation(unit, simple_idcs, complexIdx)


% Find the simple spiking neuron with the highest simple spike suppression ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ratio_base_rate_vs_suppression = zeros(1,length(simple_idcs));



% Code below selects the best match of a simple spiking neuron with the complex spiking neuron

for idx = 1:length(simple_idcs)
    
    spike_times_complex_spiking_neurons = unit.neuron(complexIdx).st;
    spike_times_simple_spiking_neurons = unit.neuron(simple_idcs(idx)).st;
    
    
    %% Finding the Cross-correlation     
    % for all the simple spiking neurons in the session with the current 
    % Complex Spiking Neuron
    
    [ST_Offsets, ~, ~] = crosscorrelogram ...
        ( spike_times_complex_spiking_neurons ...
        , spike_times_simple_spiking_neurons ...
        , getParams().cc_range ...
        );

   
   %% Finding the value of cross-correlation    
   % of every bin of the spike delay times in the cross-correlogram
    
    nST_Offsets = histcounts ...
        ( ST_Offsets ...
        , 'BinWidth', getParams().cc_BinW ...
        , 'BinLimits',getParams().cc_range ...
        , 'Normalization','probability' ...
        );
    
    %% Finding the Simple Spike Supression during the 20ms SSS window
    
    nST_Offsets_simple_spike_suppression_20ms_window = nST_Offsets(getParams().cc_minBin:getParams().cc_maxBin); % 101 to 113 are the bins corresponding to 0-20ms range of the histogram
    
    %% Finding base rate histogram
    % Histogram before the simple spike supression window to find the base rate of the histogram
    
    nST_Offsets_base_rate = nST_Offsets(1 : (getParams().cc_minBin - 1));
    
    
    %% Take out the highest (max) value
    % because sometimes the highest value is a really high strange peak and that might give an inacurate result for the cross-correlation
    
    max_ST_Offsets_idcs = sort(find(nST_Offsets == max(nST_Offsets)), 'descend');
    
    for Max_idx = max_ST_Offsets_idcs     
        if Max_idx >= getParams().cc_minBin && Max_idx <= getParams().cc_maxBin
            % If the delay time of the max value is within the 20 ms window
            % than remove it from this array
            nST_Offsets_simple_spike_suppression_20ms_window(Max_idx-(getParams().cc_minBin-1)) = [];
            
        elseif Max_idx < getParams().cc_minBin
            % If the delay time of the max value is during the base rate
            % window than remove it from this array
            nST_Offsets_base_rate(Max_idx) = [];        
            
        end
    end
    
    %% Finding the Mean of the Simple Spike Supression and the Base Rate
    mean_nST_simple_spike_supression_20ms_window = mean(nST_Offsets_simple_spike_suppression_20ms_window(1,:));
    mean_nST_base_rate = mean(nST_Offsets_base_rate);
    
    %% Finding the ratio of the Base Rate over the Simple Spike Supression
    ratio_base_rate_vs_suppression(idx) = mean_nST_base_rate...
        / mean_nST_simple_spike_supression_20ms_window;
end

%% Finding Simple Spiking Neuron with highest supression rate
% The chosen simple spike neuron is the one that gives the highest suppression rate

n_simple = simple_idcs ...
    ( find...
        ( ratio_base_rate_vs_suppression == max(ratio_base_rate_vs_suppression) ...
        , 1 ...
        ) ...
    ); 

end