function saveFigures(ax, fname)
% file = fullfile('~/Dropbox (BayesLab)/TraceExperiments/ComplexSpikeToolkit/CSpk_Suite_revised/Figures/',fname);
p = IkUtils.getParams();
file = fullfile(p.figPath,fname);
saveas(ax,file);
end