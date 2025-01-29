function loadParamTable(handles)

% paramtable = getappdata(0,'paramtable');
config = getappdata(0, 'config');

[paramfile,paramfilepath,filteridx] = uigetfile('*.csv');

if paramfile & filteridx == 1 % The filterindex thing is a hack to make sure it's a csv file
    % paramtable.data=csvread(fullfile(paramfilepath,paramfile));
    % set(handles.uitable_params,'Data',paramtable.data);
    % setappdata(0,'paramtable',paramtable);

    config.paramtable.data = csvread(fullfile(paramfilepath,paramfile));
    config.paramtable.randomize = handles.checkbox_random.Value;                 % Do we want to assume randomization here or read from gui?
    config.trialtable = makeTrialTable(config.paramtable.data, config.paramtable.randomize);

    set(handles.uitable_params, 'Data', config.paramtable.data);

    drawnow         % Seems necessary to update appdata before returning to calling function
end

setappdata(0, 'config', config);