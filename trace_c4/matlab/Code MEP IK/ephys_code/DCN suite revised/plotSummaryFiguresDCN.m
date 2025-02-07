%% I.K. 1-6-24
function [axs, figs, types] = plotSummaryFiguresDCN(mname, type, kwargs)
    arguments
        mname (1,1) string = "Shank2KOMUT"
        type (1,1) string = getMouseType(mname);
        kwargs.modulating (1,1) logical
        kwargs.event (1,1) string = "cs_facilitation"
    end


    p = IkUtils.getParams();

    switch type 
        case "Shank2KOMUT"
            sessionIdcs = [1,2]; 
            types = "Shank2KOMUT";
        case "Shank2KOWT"
            sessionIdcs = [1,2]; 
            types = "Shank2KOWT";
        otherwise
            error("Unknown session type: %s", type)
    end

    allStats = histStatsMouseWrapper ...
        ( mname ...
        , event = kwargs.event ...
        , modulating = kwargs.modulating ...
        );

    figure();

    subplot(1,2,1);
    hold on
    plot([allStats.cs.maxAmpTime],[allStats.cs.sd],'.','Color','blue');
    xlabel('Peak time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.05 0.065]);
    title_obj = title(sprintf('Peak time vs StD distribution all neurons in CS region %s',mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.cs.maxAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.cs.maxAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Peak Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryAllCS%s.png', mname), 'png')

    figure();
    non_mod_mask_cs = [allStats.cs.cs_facilitation] == 0 & [allStats.cs.cs_suppression] == 0;

    subplot(1,2,1);
    hold on
    plot([allStats.cs(non_mod_mask_cs).maxAmpTime],[allStats.cs(non_mod_mask_cs).sd],'.','Color','blue');
    xlabel('Peak time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.05 0.065]);
    title_obj = title(sprintf('Peak time vs StD distribution non-modulating CS spikes %s', mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.cs(non_mod_mask_cs).maxAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.cs(non_mod_mask_cs).maxAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Peak Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryNonModCS%s.png', mname), 'png')

    figure();
    non_mod_mask_us = [allStats.us.us_facilitation] == 0 & [allStats.us.us_suppression] == 0;

    subplot(1,2,1);
    hold on    
    plot([allStats.us(non_mod_mask_us).maxAmpTime],[allStats.us(non_mod_mask_us).sd],'.','Color','blue');
    xlabel('Peak time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.02 0.035]);
    title_obj = title(sprintf('Peak time vs StD distribution non-modulating US spikes %s', mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.us(non_mod_mask_us).maxAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.us(non_mod_mask_us).maxAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Peak Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryNonModUS%s.png', mname), 'png')

    figure();
    cs_facilitation_mask = [allStats.cs.cs_facilitation] == 1;

    subplot(1,2,1);
    hold on
    plot([allStats.cs(cs_facilitation_mask).maxAmpTime],[allStats.cs(cs_facilitation_mask).sd],'.','Color','blue');
    xlabel('Peak time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.05 0.065]);
    title_obj = title(sprintf('Peak time vs StD distribution cs facilitation spikes %s',mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.cs(cs_facilitation_mask).maxAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.cs(cs_facilitation_mask).maxAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Peak Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryCSfac%s.png', mname), 'png')

    figure();
    cs_suppression_mask = [allStats.cs.cs_suppression] == 1;

    subplot(1,2,1);
    hold on
    plot([allStats.cs(cs_suppression_mask).minAmpTime],[allStats.cs(cs_suppression_mask).sd],'.','Color','blue');
    xlabel('Dip time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.05 0.065]);
    title_obj = title(sprintf('Dip time vs StD distribution cs suppression spikes %s', mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.cs(cs_suppression_mask).minAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.cs(cs_suppression_mask).minAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Dip Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryCSsup%s.png', mname), 'png')

    figure();
    us_facilitation_mask = [allStats.us.us_facilitation] == 1;

    subplot(1,2,1);
    hold on
    plot([allStats.us(us_facilitation_mask).maxAmpTime],[allStats.us(us_facilitation_mask).sd],'.','Color','blue');
    xlabel('Peak time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.02 0.035]);
    title_obj = title(sprintf('Peak time vs StD distribution us faciliation spikes %s', mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.us(us_facilitation_mask).maxAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.us(us_facilitation_mask).maxAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Peak Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryUSfac%s.png', mname), 'png')

    figure();
    us_suppression_mask = [allStats.us.us_suppression] == 1;

    subplot(1,2,1);
    hold on
    plot([allStats.us(us_suppression_mask).minAmpTime],[allStats.us(us_suppression_mask).sd],'.','Color','blue');
    xlabel('Dip time in s')
    ylabel('Standard Deviation')
    xlim([0 0.35]); ylim([0.02 0.035]);
    title_obj = title(sprintf('Dip time vs StD distribution us suppression spikes %s', mname));
    new_position = get(title_obj, 'Position');  % Get the current position
    new_position(1) = new_position(1) + 0.15;  % Adjust the horizontal position
    set(title_obj, 'Position', new_position);  % Set the new position
    set(gca,'tickdir','out','box','off');

    subplot(1,2,2);
    hold on
    IkUtils.histogramIK([allStats.us(us_suppression_mask).minAmpTime],0:0.01:0.35,'FaceColor','blue');
    text(0.02,27,sprintf('N: %d',length([allStats.us(us_suppression_mask).minAmpTime])), "FontSize", 20);
    xlim([0 0.35]); ylim([0 30]);
    xlabel('Dip Time in s')
    set(gca,'tickdir','out','box','off');

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/SummaryUSsup%s.png', mname), 'png')

    figure()
    boxplot_data = [[allStats.cs(non_mod_mask_cs).maxAmpTime] [allStats.cs(cs_facilitation_mask).maxAmpTime] [allStats.cs(cs_suppression_mask).minAmpTime]];
    boxplot_data_cell = {[allStats.cs(non_mod_mask_cs).maxAmpTime]; [allStats.cs(cs_facilitation_mask).maxAmpTime]; [allStats.cs(cs_suppression_mask).minAmpTime]};
    grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
    boxplot_labels = {'Non modulating','CS facilitation', 'CS suppression'};
    boxplot_labels_new = {boxplot_labels{unique(grp)}};
    boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal')
    ylabel('Modulation type')
    xlabel('Peak Time')
    xlim([0 0.35]);
    set(gca,'tickdir','out','box','off');
    title(sprintf('Peak time neuron frequency per CS modulation type %s', mname))

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/boxplotCS%s.png', mname), 'png')

    figure()
    boxplot_data = [[allStats.us(non_mod_mask_us).maxAmpTime] [allStats.us(us_facilitation_mask).maxAmpTime] [allStats.us(us_suppression_mask).minAmpTime]];
    boxplot_data_cell = {[allStats.us(non_mod_mask_us).maxAmpTime]; [allStats.us(us_facilitation_mask).maxAmpTime]; [allStats.us(us_suppression_mask).minAmpTime]};
    grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
    boxplot_labels = {'Non modulating', 'US facilitation', 'US suppression'};
    boxplot_labels_new = {boxplot_labels{unique(grp)}};
    boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal')
    ylabel('Modulation type')
    xlabel('Peak Time')
    xlim([0 0.35]);
    set(gca,'tickdir','out','box','off');
    title(sprintf('Peak time neuron frequency per US modulation type %s', mname))

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/boxplotUS%s.png', mname), 'png')



    figure()
    boxplot_data = [[allStats.cs(non_mod_mask_cs).cov] [allStats.cs(cs_facilitation_mask).cov] [allStats.cs(cs_suppression_mask).cov]];
    boxplot_data_cell = {[allStats.cs(non_mod_mask_cs).cov]; [allStats.cs(cs_facilitation_mask).cov]; [allStats.cs(cs_suppression_mask).cov]};
    grp = cell2mat(arrayfun(@(i){i*ones(numel(boxplot_data_cell{i}),1)},(1:numel(boxplot_data_cell))'));
    boxplot_labels = {'Non mod','CS fac', 'CS sup'};
    boxplot_labels_new = {boxplot_labels{unique(grp)}};
    boxplot(boxplot_data, grp, 'labels', boxplot_labels_new, 'orientation','horizontal')
    ylabel('Modulation type')
    xlabel('Peak Time')
    xlim([0 0.005]);
    set(gca,'tickdir','out','box','off');
    title(sprintf('Cov neurons per CS mod type %s', mname))

    saveas(gcf, sprintf('/home/i.klinkhamer/Documents/Figures/SummaryFigures/boxplotCovCS%s.png', mname), 'png')

    
    axs = [];
    figs = [];
  

end