function ok = initializeSession()
% This is executed when user presses the button to start a new session in in the main GUI
% Initializes per session variables, some of which can be loaded from mouse config files

gui = getappdata(0, 'gui'); 
handles = guidata(gui.maingui);

metadata=getappdata(0,'metadata');
config=getappdata(0,'config');

% House keeping -- seems to be necessary to clear out persistent variables from previous sessions
clear newFrameCallbackCam1 newFrameCallbackCam2

for i=1:length(config.camera)

    metadata.cam(i).trialnum=1;
    metadata.cam(i).thresh=config.camera(i).thresh;

    metadata.cam(i).cal=0;
    metadata.cam(i).calib_offset=0;
    metadata.cam(i).calib_scale=1;
    
end

metadata.ts=[datenum(clock) 0]; % two element vector containing datenum at beginning of session and offset of current trial (in seconds) from beginning

config.userdatapath;
dirlist = dir(config.userdatapath);
dirs = [dirlist.isdir];
dirlistnames = {dirlist.name};
subdirs = dirlistnames(dirs);

mice = regexp(subdirs, '\w\w\d\d\d\d\d\d\d', 'match', 'once');
%mice = regexp(subdirs, '[a-zA-Z]\d\d\d', 'match', 'once'); % IK change. add this line from the beun versie.
mice = mice(~cellfun(@isempty, mice));

[idx, ok] = listdlg('PromptString', 'Select a mouse', 'SelectionMode', 'single', 'ListString', mice);

if ok
    mouse = mice{idx};
else
    return
end

metadata.mouse = mouse;

mousedir = fullfile(config.userdatapath, mouse); 

sessiondir = fullfile(mousedir, datestr(now,'yymmdd'));

mkdir(sessiondir)
cd(sessiondir)

metadata.basename=sprintf('%s_%s_s%02d', metadata.mouse, datestr(now,'yymmdd'), metadata.session);
metadata.folder = sessiondir;

files_in_session_folder = getFullFileNames(metadata.folder, dir(metadata.folder));

if any(contains(files_in_session_folder, metadata.basename))
    button = questdlg(sprintf('WARNING: Session directory already contains files for session number %d. \nBy continuing you will overwrite the existing files with data from the new trials.', metadata.session),...
    'WARNING - Possible data loss!', 'Continue', 'Abort Session', 'Abort Session');

    if strcmp(button, 'Abort Session')
        ok = 0;
        return
    end
end

metadata.rig = config.rig;

metadata.has_encoder = config.use_encoder;

condfile=fullfile(mousedir,'condparams.csv');

if exist(condfile, 'file')
%     paramtable = getappdata(0,'paramtable');

    config.paramtable.data = csvread(condfile);
    config.paramtable.randomize = handles.checkbox_random.Value;                 % Do we want to assume randomization here or read from gui?
    config.trialtable = makeTrialTable(config.paramtable.data, config.paramtable.randomize);

    set(handles.uitable_params, 'Data', config.paramtable.data);
    
%     setappdata(0,'paramtable',paramtable);
%     setappdata(0,'trialtable',trialtable);
end

trials.eye = [];

setappdata(0, 'trials', trials)
setappdata(0, 'metadata', metadata);
setappdata(0, 'config', config);

drawnow         % Seems necessary to update appdata before returning to calling function