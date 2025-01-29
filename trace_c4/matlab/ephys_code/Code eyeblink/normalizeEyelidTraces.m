%% IK 8-4-24
function normalized_traces = normalizeEyelidTraces(traces, behavior_trial_data)
P = IkUtils.getParams();
baseMin = P.tracesRanges.baseline.min;
baseMax = P.tracesRanges.baseline.max;
USmin = P.tracesRanges.us.min;
USmax = P.tracesRanges.us.max;
CSmin = P.tracesRanges.cs.min;

baseline_min = min(nanmean(traces(:, baseMin:baseMax), 2));      % Calculate baselines as the mean of traces within the specified range

meanUR = nanmean(max(traces(:,USmin:USmax)'));

fullClosureRange = meanUR - baseline_min;

maskEyeMoreThanHalfOpen = nanmean(traces(:,baseMin:baseMax)') < 0.5 * fullClosureRange;

traces(~maskEyeMoreThanHalfOpen,:) = nan;

traces_2baseline = traces - ...
    (nanmean(traces(:,baseMin:baseMax), 2) - ...
    baseline_min);

baseline_new_min = min(...
    nanmean(traces_2baseline(behavior_trial_data.c_csdur==0,baseMin:baseMax)')); %#ok<*UDIM> 

stdev = nanstd(traces_2baseline(:,baseMin:baseMax), 0, 2); %#ok<*NANSTD> 

mask_stdev = stdev < 3*nanmean(stdev); %#ok<*NANMEAN> 
mask_lowest = min(traces_2baseline(:,CSmin:end), [], 2) > (baseline_new_min - 0.1);
mask_highest = max(traces_2baseline(:,USmax+1:end), [],2) < max(max(traces_2baseline(:,USmin:USmax)'));
traces_2baseline(~mask_lowest | ~mask_highest|~mask_stdev, :) = nan;

min_val = min(nanmean(traces_2baseline(:,baseMin:baseMax)'));
max_val = max(max(traces_2baseline(:,USmin:USmax)'));
normalized_traces = (traces_2baseline - min_val) / (max_val - min_val);
end