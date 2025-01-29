%% IK 15-3-2024     Overview eyeblink training and eyelid traces over the course of training
function overviewTraces()
clear
close all
app = getUImain(field = "ProcessonlyEphysmiceCheckBox", value = false, scriptname=mfilename);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        % If the UI figure has been closed, exit the loop
        break;
    end
    disp("Calculating...")
    appData = evalin('base', 'appData');

    showBGlines = appData.showBGlines;
    savefigures = appData.savefigures; % Set this to false if it's taking too long to plot.
    day = appData.behaviorLastDay;
    Gene = appData.gene;
    loc = appData.loc;

    if appData.enumeration == "All mice of one gene" && appData.gene == "Shank2"
        disp("Not enough Shank2 L7 mice to do this. Please pick 'All mice of both types of a gene + location' with 'Shank2' and 'KO' and try again.")
        continue
    end

    mouseCodes = appData.mousecodes;
    mouseTypes = appData.mousetypes;
    non_valid_idcs = find(mouseCodes == "21-MI10159-01"|mouseCodes == "21-MI10159-02"|mouseCodes == "21-MI10159-06"|mouseCodes == "Shank2KOMUT"|mouseCodes == "Shank2KOWT");
    mouseCodes(non_valid_idcs) = [];
    mouseTypes(non_valid_idcs) = [];
    if isempty(mouseCodes)
        disp("No eyelid data found. Pick another mouse.")
        continue
    end

    if appData.enumeration == "All mice" || appData.enumeration == "Single mouse"
        for m = 1:length(mouseCodes)
            close all
            dataStruct = calcMouseForOverviewTraces(mouseCodes(m)); % calculations
            plotEyeTracesFiguresSingleMouse(dataStruct, mouseCodes(m), savefigures, day) % Final figures
            plotCRoverviewSingleMouse(dataStruct, mouseCodes(m), savefigures, day)
        end
    elseif appData.enumeration == "All mice of one gene" || appData.enumeration == "All mice of both types of a gene + location"
        dataStruct = struct;
        selectedMiceTypes = sort(unique(appData.mousetypes));
        types = IkUtils.getParams().types;

        for t = 1:length(selectedMiceTypes)
            t_parts = split(selectedMiceTypes(t), ' ');
            mouseTypeCC = strjoin(t_parts, '');
            idx = find([types.full] == mouseTypeCC);
            Type = types.type(idx);
            Gene = types.gene(idx);
            Loc = types.loc(idx);
            dataStruct(t).type = selectedMiceTypes(t);
            mouseListType = mouseCodes(mouseTypes == selectedMiceTypes(t));
            dataStruct(t).data = calcMouseForOverviewTraces(mouseListType, Gene, Loc, Type);
        end
       
        %% Final figures
        dataStructTypes = [dataStruct.type];
        if appData.enumeration == "All mice of both types of a gene + location"            
            WT = dataStruct(contains(dataStructTypes, "WT")).data;
            MUT= dataStruct(contains(dataStructTypes, "MUT")).data;

            plotEyeTracesFigures(WT, MUT, Gene, loc, savefigures, day)
            plotCRoverview(WT, MUT, Gene, loc, savefigures, showBGlines, day)
        else
            WTKO = dataStruct(contains(dataStructTypes, "WT") & contains(dataStructTypes, "KO")).data;
            MUTKO = dataStruct(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "KO")).data;
            WTL7 = dataStruct(contains(dataStructTypes, "WT") & contains(dataStructTypes, "L7")).data;
            MUTL7 = dataStruct(contains(dataStructTypes, "MUT") & contains(dataStructTypes, "L7")).data;

            plotEyeTracesFigures(WTKO, MUTKO, Gene, "KO", savefigures, day)
            plotEyeTracesFigures(WTL7, MUTL7, Gene, "L7", savefigures, day)
            plotCRoverviewAllGroups(WTKO, MUTKO, WTL7, MUTL7, Gene, savefigures, showBGlines, day)
            plotOnsetLatencyAllGroups(WTKO, MUTKO, WTL7, MUTL7, Gene, savefigures, day)
        end
    end

    disp("End of loop. Waiting for new UI input.")
end
end

function plotEyeTracesFiguresSingleMouse(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotEyeTracesFiguresSingleMouse')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end

eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotEyelidTraces(Mouse, mcode, eyeblinkTracesOverview, [2,1,1], Day)
plotEyeTracesOverTime(Mouse, mcode, eyeblinkTracesOverview, [2,1,2], Day)

if savefigures
    fname = sprintf('eyeblinkTracesOverview_%s.eps', mcode);
    file = fullfile(IkUtils.getParams().figPath, "Training overview eps",fname);
%     saveas(eyeblinkTracesOverview,file);
    fig2save = eyeblinkTracesOverview;
    print(fig2save, '-depsc', '-painters', file)

    fname2 = sprintf('eyeblinkTracesOverview_%s.png', mcode);
    file = fullfile(IkUtils.getParams().figPath, "Training overview",fname2);
    saveas(eyeblinkTracesOverview,file);
end
end

function plotCRoverviewSingleMouse(Mouse, mcode, savefigures, Day)
P = IkUtils.getParams();
if nargin < nargin('plotCRoverviewSingleMouse')
    Day = min(IkUtils.getParams().lastDay, size(Mouse.onsets,1));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure(Mouse, mcode, CROverview, [2,2,1], day = Day);
plotCRperc(Mouse, mcode, CROverview, [2,2,2], day = Day);
plotOffset(Mouse, mcode, CROverview, [2,2,3], Day);

if savefigures
    fname = sprintf('CRoverview_%s.eps', mcode);
    file = fullfile(IkUtils.getParams().figPath, "Training overview eps",fname);
%     saveas(CROverview,file);
        fig2save = CROverview;
    print(fig2save, '-depsc', '-painters', file)
    fname2 = sprintf('CRoverview_%s.png', mcode);
    file = fullfile(IkUtils.getParams().figPath, "Training overview",fname2);
    saveas(CROverview,file);
end
end

function plotEyeTracesFigures(WT, MUT, Gene, loc, savefigures, Day) %#ok<*DEFNU> 
if nargin < nargin('plotEyeTracesFigures')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
eyeblinkTracesOverview = figure('Color', 'white','Position', [10 0 1000 700]);

WT_CRamp = squeeze(WT.cramp(Day, :, :));
WT_CRamp_flat = WT_CRamp(:);
WT_CRamp_mean = nanmean(WT_CRamp, 2);
WT_CRamp_mean = WT_CRamp_mean(~isnan(WT_CRamp_mean));
MUT_CRamp = squeeze(MUT.cramp(Day, :, :));
MUT_CRamp_flat = MUT_CRamp(:);
MUT_CRamp_mean = nanmean(MUT_CRamp, 2);
MUT_CRamp_mean = MUT_CRamp_mean(~isnan(MUT_CRamp_mean));

[h, p, chi, stats] = ttest2(WT_CRamp_flat, MUT_CRamp_flat)
[h_mean, p_mean, chi_mean, stats_mean] = ttest2(WT_CRamp_mean, MUT_CRamp_mean)

plotEyelidTraces(WT, [Gene, loc, "WT"], eyeblinkTracesOverview, [2,2,1], Day)
plotEyelidTraces(MUT, [Gene, loc, "MUT"], eyeblinkTracesOverview, [2,2,2], Day)
plotEyeTracesOverTime(WT, [Gene, loc, "WT"], eyeblinkTracesOverview, [2,2,3], Day)
plotEyeTracesOverTime(MUT, [Gene, loc, "MUT"], eyeblinkTracesOverview, [2,2,4], Day)

if savefigures
    fname = sprintf('eyeblinkTracesOverview%s%s_%s.eps', Gene, loc, IkUtils.Now);
%     saveFigures(eyeblinkTracesOverview, fname);
    fig2save = eyeblinkTracesOverview;
    print(fig2save, '-depsc', '-painters', fullfile(IkUtils.getParams().figPath, 'Training overview groups', fname))
    fname2 = sprintf('eyeblinkTracesOverview%s%s_%s.png', Gene, loc, IkUtils.Now);
    saveFigures(eyeblinkTracesOverview, fname2);
end
end

function plotCRoverviewAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, showBGlines, Day)
if nargin < nargin('plotCRoverviewAllGroups')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotErrorBarFigure([WT,MUT], [Gene, "KO"], CROverview, [2,2,1], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WT, MUT], [Gene, "KO"], CROverview, [2,2,2], showBackGroundLines = showBGlines, day = Day);
plotErrorBarFigure([WTL7,MUTL7], [Gene, "L7"], CROverview, [2,2,3], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WTL7, MUTL7], [Gene, "L7"], CROverview, [2,2,4], showBackGroundLines = showBGlines, day = Day);

if savefigures
    fname = sprintf('CRpercentage%s_%s.eps', Gene, IkUtils.Now);
%     saveFigures(CROverview, fname);
        fig2save = CROverview;
    print(fig2save, '-depsc', '-painters', fullfile(IkUtils.getParams().figPath, 'Training overview groups', fname))
    fname2 = sprintf('CRpercentage%s_%s.png', Gene, IkUtils.Now);
    saveFigures(CROverview, fname2);
end
end

function plotOnsetLatencyAllGroups(WT, MUT, WTL7, MUTL7, Gene, savefigures, Day)
if nargin < nargin('plotOnsetLatencyAllGroups')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
LatencyOverview = figure('Color', 'white','Position', [10 0 1000 700]);

plotOffset(WT, [Gene, "KO", "WT"], LatencyOverview, [2,2,1], Day)
plotOffset(MUT, [Gene, "KO", "MUT"], LatencyOverview, [2,2,2], Day)
plotOffset(WTL7, [Gene, "L7", "WT"], LatencyOverview, [2,2,3], Day)
plotOffset(MUTL7, [Gene, "L7", "MUT"], LatencyOverview, [2,2,4], Day)

if savefigures
    fname = sprintf('CRonsetLatency%s_%s.eps', Gene, IkUtils.Now);
%     saveFigures(LatencyOverview, fname);
            fig2save = LatencyOverview;
    print(fig2save, '-depsc', '-painters', fullfile(IkUtils.getParams().figPath, 'Training overview groups', fname))
    fname2 = sprintf('CRonsetLatency%s_%s.png', Gene, IkUtils.Now);
    saveFigures(LatencyOverview, fname2);
end
end

function plotCRoverview(WT, MUT, Gene, loc, savefigures, showBGlines, Day)
if nargin < nargin('plotCRoverview')
    Day = min(IkUtils.getParams().lastDay, min(size(WT.onsets,1),size(MUT.onsets,1)));
end
CROverview = figure('Color', 'white','Position', [10 0 1000 700]);

onsets_firstday_WT_ = WT.onsets(1, :);% Extract onset times for WT and MUT
onsets_firstday_MUT_ = MUT.onsets(1, :);% Extract onset times for WT and MUT
onsets_lastday_WT_ = WT.onsets(Day, :);
onsets_lastday_MUT_ = MUT.onsets(Day, :);

onsets_firstday_WT = onsets_firstday_WT_(~isnan(onsets_firstday_WT_));% Remove NaN values
onsets_lastday_WT = onsets_lastday_WT_(~isnan(onsets_lastday_WT_));
onsets_firstday_MUT = onsets_firstday_MUT_(~isnan(onsets_firstday_MUT_));% Remove NaN values
onsets_lastday_MUT = onsets_lastday_MUT_(~isnan(onsets_lastday_MUT_));

[h, p] = kstest2(onsets_firstday_WT, onsets_firstday_MUT);% Perform KS test
disp(['KS test result first day WT vs. MUT: h = ', num2str(h), ', p = ', num2str(p)])% Display KS test result
[h, p] = kstest2(onsets_lastday_WT, onsets_lastday_MUT);% Perform KS test
disp(['KS test result last day WT vs. MUT: h = ', num2str(h), ', p = ', num2str(p)])% Display KS test result

% Perform Mann-Whitney U test
[p_mw, h_mw] = ranksum(onsets_firstday_WT, onsets_firstday_MUT);
disp(['Mann-Whitney U test result first day WT vs. MUT: h = ', num2str(h_mw), ', p = ', num2str(p_mw)])
% Perform Mann-Whitney U test
[p_mw, h_mw] = ranksum(onsets_lastday_WT, onsets_lastday_MUT);
disp(['Mann-Whitney U test result last day WT vs. MUT: h = ', num2str(h_mw), ', p = ', num2str(p_mw)])

% % Perform Wilcoxon Signed-Rank Test
% onsets_firstday_WT = onsets_firstday_WT_(~isnan(onsets_firstday_WT_) & ~isnan(onsets_lastday_WT_));
% onsets_lastday_WT = onsets_lastday_WT_(~isnan(onsets_firstday_WT_) & ~isnan(onsets_lastday_WT_));
% onsets_firstday_MUT = onsets_firstday_MUT_(~isnan(onsets_firstday_MUT_) & ~isnan(onsets_lastday_MUT_));
% onsets_lastday_MUT = onsets_lastday_MUT_(~isnan(onsets_firstday_MUT_) & ~isnan(onsets_lastday_MUT_));
% [p_ws, h_ws] = signrank(onsets_firstday_WT, onsets_firstday_MUT);
% disp(['Wilcoxon Signed-Rank Test First Day WT vs. MUT: h = ', num2str(h_ws), ', p = ', num2str(p_ws)]);
% [p_ws, h_ws] = signrank(onsets_lastday_WT, onsets_lastday_MUT);
% disp(['Wilcoxon Signed-Rank Test Last Day WT vs. MUT: h = ', num2str(h_ws), ', p = ', num2str(p_ws)]);

plotErrorBarFigure([WT,MUT], [Gene, loc], CROverview, [2,2,1], showBackGroundLines = showBGlines, day = Day);
plotCRperc([WT, MUT], [Gene, loc], CROverview, [2,2,2], showBackGroundLines = showBGlines, day = Day);
plotOffset(WT, [Gene, loc, "WT"], CROverview, [2,2,3], Day);
plotOffset(MUT, [Gene, loc, "MUT"], CROverview, [2,2,4], Day);

if savefigures
    fname = sprintf('CRoverview%s%s_%s.eps', Gene, loc, IkUtils.Now);
%     saveFigures(CROverview, fname);
          fig2save = CROverview;
    print(fig2save, '-depsc', '-painters', fullfile(IkUtils.getParams().figPath, 'Training overview groups', fname))
    fname2 = sprintf('CRoverview%s%s_%s.png', Gene, loc, IkUtils.Now);
    saveFigures(CROverview, fname2);
end
end