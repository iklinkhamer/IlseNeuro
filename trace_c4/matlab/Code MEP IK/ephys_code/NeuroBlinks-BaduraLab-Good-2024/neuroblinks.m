% Set up environment and launch app 
function neuroblinks(varargin)
    
    [basedir, mfile, ext] = fileparts(mfilename('fullpath'));
    oldpath = addpath(genpath_exclude(fullfile(basedir), '.git'));

    % Load local configuration for these rigs
    % Should be somewhere in path but not "neuroblinks" directory or subdirectory
    neuroblinks_config
    rig_config;
    user_config;

    % Set up defaults in case user doesn't specify all options
    rig = DEFAULTRIG;
    device = DEFAULTDEVICE;

    % Special case if we just want to set up path and configuration but not launch full program
    if length(varargin) > 0 && varargin{1} == 0
        setappdata(0, 'config', config)
        return
    end

    % Input parsing
    % Values passed override defaults set in neuroblinks_config
    if nargin > 0
        for i = 1:nargin
            if any(strcmpi(varargin{i}, config.ALLOWEDDEVICES))
                device = varargin{i};
            elseif ismember(varargin{i}, 1:length(config.camera(1).IDS))
                rig = varargin{i}; 
            end
        end
    end

    config.rig = rig;
    config.device = device;
    config.oldpath = oldpath;   % save so we can restore later if desired

    setappdata(0, 'config', config)

    drawnow         % Seems necessary to update appdata before calling new function
    
    launch(config)