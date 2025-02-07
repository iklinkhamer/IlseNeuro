function processSessions(folder)
% Traverse directories starting at ROOTFOLDER and recursively check for compressed video subdirectory, running processTrials when found.

TEST = 0;

% Check current directory for subdirectory containing the compressed files and if there is none, recursively go to the next nested directory
dirInfo= dir(folder);

isDir=[dirInfo.isdir];
dirNames={dirInfo(isDir).name};
dirNames(strcmp(dirNames, '.') | strcmp(dirNames, '..'))=[];

if isempty(dirNames)
	return	% This will allow us to stop recursion
end

for i=1:length(dirNames)
	% If the directory has 6 digits assume this is a session directory
	if containsregexp(dirNames{i},'\d\d\d\d\d\d')
		disp(['Found folder ' fullfile(folder, dirNames{i})])
		if ~TEST
			% This is where we do what it is we want to do
			% Uncomment the specific lines you want to run on each directory

			% Run process trials on raw videos
			if ~exist(fullfile(folder, dirNames{i}, 'trialdata.mat'),'file')
				trials = processTrials(fullfile(folder, dirNames{i}), 'recalibrate');	% Run processTrials on base session directory
				save(fullfile(folder, dirNames{i}, 'trialdata.mat'), 'trials')
			end

			% % Run process trials on compressed videos
			% if ~exist(fullfile(folder, dirNames{i}, 'trialdata.mat'),'file')
			% 	trials = processTrials(fullfile(folder, dirNames{i}, 'compressed'), 'recalibrate');	% Run processTrials on base session directory
			% 	save(fullfile(folder, dirNames{i}, 'trialdata.mat'), 'trials')
			% end

			% % Put encoder data in trial structure
			% if exist(fullfile(folder, dirNames{i}, 'trialdata.mat'),'file')
			% 	t = load(fullfile(folder, dirNames{i}, 'trialdata.mat'));
			% 	trials = addEncoderData(t.trials, fullfile(folder, dirNames{i}));
			% 	save(fullfile(folder, dirNames{i}, 'trialdata.mat'), 'trials')
			% 	fprintf('Added encoder data to %s\n', fullfile(folder, dirNames{i}, 'trialdata.mat'))
			% end

		end
	else
		processSessions(fullfile(folder,dirNames{i}))
	end

end

function bool = containsregexp(str, pattern)

	bool = ~isempty(regexp(str, pattern, 'match', 'once'));