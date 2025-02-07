function [filtered, filterMask] = rejectSDOutliers(dataPoints, kwargs)
    arguments
        dataPoints
        kwargs.plotResult (1,1) logical = false;
    end
    
    center = mean(dataPoints);
    sd = std(dataPoints);
    
    centerDistances = abs(dataPoints - center);
    filterMask = centerDistances <= 3*sd;
    
    filtered = dataPoints(filterMask);
    
    if kwargs.plotResult
        ax = IkUtils.initPlots([1 1]);
    
        plot(ax, dataPoints, 1, '*r')
        plot(ax, filtered, 1, '*b')
    end
    
    
end