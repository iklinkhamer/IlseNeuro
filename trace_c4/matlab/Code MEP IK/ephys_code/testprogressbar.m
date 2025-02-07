fprintf("|----"); % Print initial bar
prevPercent = 0;
maxIt = 100;

for i = 1:maxIt
    % add your loop here

    nTotalDataSamples = maxIt;
    percentComplete = i / nTotalDataSamples * 100;
    if floor(percentComplete) > prevPercent || percentComplete == 100
        fprintf(repmat('\b', 1, numel(num2str(floor(percentComplete))) + 3));
        fprintf('*');
        fprintf('| %d%%', floor(percentComplete));
        prevPercent = floor(percentComplete) + 1;
    end

    pause(0.05);
end
