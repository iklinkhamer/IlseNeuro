function rasterFiltered = rasterFilter(raster, range)
if nargin < 2
    range.min = min(raster(1,:));
    range.max = max(raster(1,:));
end
binWidth = IkUtils.getParams().BinW;
nTrials = max(raster(2,:)) - IkUtils.getParams().raster_spike_height;
binEdges = range.min:binWidth:range.max;
filterMask = zeros(1,length(raster));
for t = 1 : nTrials
    spikeTimesTrial_ = raster(1,raster(2,:) == t);
    rangeMask = ...
        spikeTimesTrial_ >= range.min ...
        & spikeTimesTrial_ <= range.max;
    spikeTimesTrial = spikeTimesTrial_(rangeMask);
    amplitudes = histcounts(spikeTimesTrial, binEdges);
    baseline = mean(amplitudes);
    %rate = mean(amplitudes/nTrials/binWidth);
    baserate = baseline/binWidth;
    % nSpikesTrial(t) = length(find(raster(2,:) == t));
    if baserate > IkUtils.getParams().Sspk_base_filter %nSpikesTrial(t) > 0
        filterMask(raster(2,:) == t) = 1;
        filterMask(find(raster(2,:) == t) + 1) = 1;
        filterMask(find(raster(2,:) == t) + 2) = 1;
%     else
%         fprintf("Trial %d baserate too low: %d\n", t, baserate);
    end
end
rasterFiltered(1,:) = raster(1,find(filterMask));
rasterFiltered(2,:) = raster(2,find(filterMask));
rasterFiltered(3,:) = NaN;
rasterFiltered(3, rasterFiltered(2,:) <= 20.9) = 0;
rasterFiltered(3, rasterFiltered(2,:) > 20.9 & rasterFiltered(2,:) <= 220.9) = 1;
rasterFiltered(3, rasterFiltered(2,:) > 220.9) = 2;
end