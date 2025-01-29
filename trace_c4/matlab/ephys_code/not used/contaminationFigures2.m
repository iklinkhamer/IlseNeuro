function contaminationFigures2(ax, contamination)
    
    histogram ...
        ( ax ...
        , contamination.spike_intervals ...
        , 'BinWidth', getParams().cBinW ...
        , 'BinLimits', getParams().cBinLimits...
        , 'EdgeColor', 'none' ...
        , 'FaceColor', getColors().c_blue ...
        );

    % Contamination is Difference between spike times consecutive spikes
    t = title(ax, sprintf('Contamination = %.2f',contamination.percentage));
    if contamination.percentage > getParams().cThreshold
        set(t, 'Color', "r");
    else
        set(t, 'Color', "k");
    end

    box(ax, 'off');
    set(ax, 'TickDir','out');
    xlim(ax, [0 1200]);
    % ylim([0 100]);
    xlabel(ax, 'Difference in spike time in ms');
    ylabel(ax, 'Amount of contamination');


end