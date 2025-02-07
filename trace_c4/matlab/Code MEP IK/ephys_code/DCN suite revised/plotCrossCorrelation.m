% TODO: This function has thrown an error during curation, find the bug and/or
% wrap its calling in a try/catch block
function plotCrossCorrelation(f, unit, bestSimpleIdx, complexIdx, mname, s) 
%%  Cross-correlation figures

try
clf(f)
catch
end

if isempty(bestSimpleIdx)
    return;
end

figure(f);
%axs1 = IkUtils.initPlots([1,1], parents = fig(1));

spike_times_n_simple = unit.neuron(bestSimpleIdx).st;
spike_times_complex_spiking_neurons = unit.neuron(complexIdx).st;

[ST_Offsets_n_simple, ~, ~] = crosscorrelogram...
    (spike_times_complex_spiking_neurons...
    ,spike_times_n_simple...
    ,getParams().cc_range...
    );

histogram...
    ( ST_Offsets_n_simple...
    ,'BinWidth', getParams().cc_BinW...
    ,'BinLimits',getParams().cc_range...
    ,'Normalization','probability'...
    ,'EdgeColor','none'...
    ,'FaceColor',getColors().c_blue...
    );

hold on

title...
    (sprintf('Cross-correlation %s', mname)...
    ,sprintf('S: %d, N: %d with N simple: %d'...
    ,s...
    ,unit.KSTrueID(complexIdx)...
    ,unit.KSTrueID(bestSimpleIdx))...
    );

xlim(getParams().cc_range);
xticks([-0.2 -0.1 0 0.1 0.2]);
xticklabels([-200 -100 0 100 200]);

xlabel('Lags in ms');
ylabel('Amount of correlation');

% Figure settings
set(0, 'defaultFigureRenderer', 'painters');
set(gca,'TickDir','out');
box('off');
box(gca,'off');

hold off
end