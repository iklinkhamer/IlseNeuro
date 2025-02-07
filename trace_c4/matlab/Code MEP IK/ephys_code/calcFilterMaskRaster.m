function rasterFiltered = calcFilterMaskRaster(raster)
nTrials = floor(max(raster(2,:)));
filterMask = zeros(1,length(raster));
for t = 1 : nTrials
    nSpikesTrial(t) = length(find(raster(2,:) == t));
    if nSpikesTrial(t) > 0
        filterMask(raster(2,:) == t) = 1;
        filterMask(find(raster(2,:) == t) + 1) = 1;
        filterMask(find(raster(2,:) == t) + 2) = 1;
    end
end

rasterFiltered(1,:) = raster(1,find(filterMask));
rasterFiltered(2,:) = raster(2,find(filterMask));

end