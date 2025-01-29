base_dir = 'C:\Users\nmhet\Documents\studie\MEP\data';

%%
tic
mouse = 'MI1544101';
isi = 250;
session = 2;
us = 3;
cs = 1;
% day = '210308';

day_offset = 1;
day = datestr(now - day_offset, 'yymmdd')

folder = fullfile(base_dir,  mouse, day);

trials = processTrials(folder) %,'recalibrate');  % Recalibrate eyelid

%%
if ~isempty(trials)
    
    save(fullfile(folder, 'trialdata.mat'), 'trials');

    baseline = mean(mean(trials.eyelidpos(trials.c_csdur==0,1:40)));
    fullclosure = mean(max(trials.eyelidpos(trials.c_csdur==0,40:95),[],2));
    
    traces = (trials.eyelidpos - baseline) / (fullclosure - baseline);
    
    figure 
    plot(trials.tm(1,:), traces)
%     hf1 = makePlots(trials, isi, session, us, cs);
% 
%     hgsave(hf1, fullfile(folder, 'CRs.fig'));
%    % hgsave(hf2, fullfile(folder, 'CR_amp_trend.fig'));
% 
%     print(hf1, fullfile(folder, sprintf('%s_%s_CRs.pdf', mouse, day)), '-dpdf')
%    % print(hf2, fullfile(folder, sprintf('%s_%s_CR_amp_trend.pdf', mouse, day)), '-dpdf')

end

toc