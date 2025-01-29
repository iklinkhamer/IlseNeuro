p = IkUtils.getParams();
timestamps = p.psthRanges.cs_full.min:p.spikeTimesBinW:p.psthRanges.cs_full.max;
stimlabels = zeros(length(timestamps),4);
stimlabels(:,1) = timestamps < p.sspkRanges.cs.min;
stimlabels(:,2) = timestamps >= p.sspkRanges.cs.min & timestamps < p.sspkRanges.cs.max;
stimlabels(:,3) = timestamps >= p.sspkRanges.us.min & timestamps < p.sspkRanges.us.max;
stimlabels(:,4) = timestamps >= p.sspkRanges.us.max;
save(fullfile("/home/i.klinkhamer/Documents/Data/",'stimlabels.mat'),'stimlabels')