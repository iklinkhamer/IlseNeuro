% Compute contamination for a single data channel
function contamination = computeContamination(channelData)
    
    spikeTimes = channelData.st;

    contamination.spike_intervals = diff(spikeTimes) * 1000;

    shortIntervals = sum(contamination.spike_intervals < 5);
    
    totalIntervals = numel(contamination.spike_intervals);
    
    contamination.percentage = shortIntervals / totalIntervals * 100;
    
end