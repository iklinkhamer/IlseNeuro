function trialvars=readTrialTable(current_tr)
% Return current row from trial table, which is stored in the app data for the root figure

config = getappdata(0, 'config');

[m,n] = size(config.trialtable);

current_tr=mod(current_tr-1,m)+1;	% Cycle through the table again if we've reached the end

trialvars=config.trialtable(current_tr,:);
