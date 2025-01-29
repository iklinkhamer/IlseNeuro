function trials=processTrials(folder, varargin)
% TRIALS=processConditioningTrials(FOLDER)
% Return trials structure containing eyelid data and trial parameters for all trials in a session
% Optionally, specify threshold for binary image and maximum number of video frames per trial to use for extracting eyelid trace
% Note that this is a simplification of prior versions that required calibration info as an argument so it's not compatible with 
% function calls that assume additional arguments. This version should be more universally useful because it uses the calibration
% info stored in the metadata.


% if length(varargin) > 0
% 	thresh=varargin{1};
% end

% % This part is sloppy but I don't want to change all the other logic right now
% USESAVEDCALIB = 0;

% % Error checking
% if isstruct(calib)
% 	if ~isfield(calib,'scale') || ~isfield(calib,'offset')
% 		error('You must specify a valid calibration structure or file from which the structure can be computed')
% 	end
% elseif exist(calib,'file')
% 	[data,metadata]=loadCompressed(calib);
% 	if ~exist('thresh')
% 		thresh=metadata.cam.thresh;
% 	end
% 	[y,t]=vid2eyetrace(data,metadata,thresh,5);
% 	calib=getcalib(y);
% elseif isempty(calib)
% 	% use the calibration stored in each file's metadata
% 	USESAVEDCALIB = 1;
% else
% 	error('You must specify a valid calibration structure or file from which the structure can be computed')
% end

% By now we should have a valid calib structure to use for calibrating all files

ms2frm = @(ms) ms ./ 1e3 * 200;	% Convert ms to frames

if length(varargin) > 0 	% If user specifies 'recalibrate as argument', recalbrate the eyelid using global min and max during trial (up to time of UR to exclude paw artifacts) 
	if strcmp(varargin{1}, 'recalibrate')
		RECALIBRATE = 1;
	else
		RECALIBRATE = 0;
    end
else
    RECALIBRATE = 0;
end

if ~exist(folder,'dir')
	error('The directory you specified (%s) does not exist',folder);
end

% Get our directory listing
% fnames=getFullFileNames(folder,dir(fullfile(folder,'*.avi')));		% for compressed AVI files
% fnames = getFullFileNames(folder, dir(fullfile(folder,'*.mp4')));		% for compressed MP4 files
fnames = getFullFileNames(folder, dir(fullfile(folder,'*.mat')));		% for compressed MAT files

% Only keep the files that match the pattern MOUSEID_DATE_SXX or MOUSEID_DATE_TXX, skipping for instance trials from Camera 2 or encoder
matches = regexp(fnames, '\w\w\d\d\d\d\d\d\d_\d\d\d\d\d\d_[ts]\d\d[a-z]?_\d\d\d_cam1', 'start', 'once');
fnames = fnames(cellfun(@(a) ~isempty(a), matches));

if isempty(fnames)
	warning('No files found in directory %s', folder);
	trials = [];
	return
end

% Preallocate variables so we can use parfor loop to process the files
eyelidpos = cell(length(fnames), 1);	% We have to use a cell array because trials may have different lengths
tm = cell(length(fnames), 1);			% Same for time

% encoder_time_cell = cell(length(fnames), 1);			% Assume encoder has same time base as eyelid so don't need to save time
encoder_displacement_cell = cell(length(fnames), 1);	% We have to use a cell array because trials may have different lengths

c_isi = NaN(length(fnames), 1);
c_csnum = NaN(length(fnames), 1);
c_csdur = NaN(length(fnames), 1);
c_csintensity = NaN(length(fnames), 1);
c_usdur = NaN(length(fnames), 1);
c_usnum = NaN(length(fnames), 1);
c_whitenoise = NaN(length(fnames), 1);

trialnum = zeros(length(fnames), 1);
ttype = cell(length(fnames), 1);

trialtime = zeros(length(fnames), 1);

numframes = zeros(length(fnames), 1);

% Use a parallel for loop to speed things up
% matlabpool open	% Start a parallel computing pool using default number of labs (usually 4-8).
% cleaner = onCleanup(@() matlabpool('close'));

parfor i=1:length(fnames)

	[p,fname,ext] = fileparts(fnames{i});

	basename = regexp(fname, '\w\w\d\d\d\d\d\d\d_\d\d\d\d\d\d_[ts]\d\d[a-z]?_\d\d\d', 'match', 'once');

	try	
		% We now prefilter above so the only condition that can be true is based on how we filtered extensions above
		% Leaving this in for flexibility if we decide to change the filter or something
		% Have to prefix everything with d.* so that load will work in parfor
		if contains(ext, [".mp4", ".avi"])			
			[d.vid,d.metadata] = loadCompressed(fnames{i});
		else
			d = load(fnames{i});
		end
    catch
        fprintf('Problem with file %s', fnames{i})
	end

    calib = struct;
 
	calib.scale = d.metadata.cam(1).calib_scale;
	calib.offset = d.metadata.cam(1).calib_offset;

	thresh = d.metadata.cam(1).thresh; 
	
	if RECALIBRATE	% if we're going to recalibrate (at end of function) just get raw pixel counts here
		[eyelidpos{i},tm{i}] = vid2eyetrace(d.vid, d.metadata, thresh, 5);
	else
		[eyelidpos{i},tm{i}] = vid2eyetrace(d.vid, d.metadata, thresh, 5, calib);
	end
	
	tm{i} = tm{i} .* 1e3; 		% convert to ms

	c_isi(i) = d.metadata.stim.c.isi;
	c_csnum(i) = d.metadata.stim.c.csnum;
	c_csdur(i) = d.metadata.stim.c.csdur;
	c_csintensity(i) = d.metadata.stim.c.csintensity;
	c_usnum(i) = d.metadata.stim.c.usnum;
	c_usdur(i) = d.metadata.stim.c.usdur;
    
    if isfield(d.metadata.stim.c, 'whitenoise')
        c_whitenoise(i) = d.metadata.stim.c.whitenoise;
    end

	encoder_fname = fullfile(p, sprintf('%s_encoder.mat', basename));

	if exist(encoder_fname, 'file')

		e = load(encoder_fname);

		% encoder_time_cell{i} = e.encoder.time;
		encoder_displacement_cell{i} = e.encoder.displacement;

	end
	
	trialnum(i) = d.metadata.cam(1).trialnum;
	ttype{i} = d.metadata.stim.type;

	trialtime(i) = d.metadata.ts(2);

	numframes(i) = length(eyelidpos{i});

	fprintf('Processed file %s\n', fname)
	
end

disp('Done reading data')

% matlabpool close


MAXFRAMES = max(numframes);

% Now that we know how long each trial is turn the cell arrays into matrices

eyelid_traces = NaN(length(fnames), MAXFRAMES);
eyelid_times = NaN(length(fnames), MAXFRAMES);

% We assume that encoder has same sampling rate as eyelid
encoder_displacement = NaN(length(fnames), MAXFRAMES);
% encoder_time = NaN(length(fnames), MAXFRAMES);

try
	for i=1:length(fnames) 

		eyelid_trace = eyelidpos{i}; 
		t = tm{i}; 

		enc_displacement = encoder_displacement_cell{i};
		% enc_time = encoder_time_cell{i};

		en = length(eyelid_trace); 
		if en > MAXFRAMES
			en = MAXFRAMES; 
		end 

		eyelid_traces(i,1:en) = eyelid_trace(1:en); 
		eyelid_times(i,1:en) = t(1:en); 

		encoder_displacement(i,1:en) = enc_displacement(1:en);
		% encoder_time(i,1:en) = enc_time(1:en);
		
	end
catch
    disp(i)
end

if RECALIBRATE		% Recalbrate eyelid traces by taking global min and max of traces during all trials (within expected response period)
    % Note that this will only work if there are US or paired trials in
    % session
	pretm = 200;
	pre = 1:ms2frm(pretm);
	win = (ms2frm(pretm) + ms2frm(mode(c_isi)):ms2frm(pretm) + ms2frm(mode(c_isi) + 50)) + 1;
	min_el = min(mean(eyelid_traces(c_usdur > 0,pre), 1));
	max_el = max(mean(eyelid_traces(c_usdur > 0,win), 1));

	% y=(tr-calib.offset)./calib.scale;
	eyelid_traces = (eyelid_traces - min_el) ./ (max_el - min_el);
end

session_cell = regexp(fnames, '_s(\d\d)[a-d]?_', 'tokens', 'once');
session_number = cellfun(@safestr2double, session_cell);

trials.session_of_day = session_number;

trials.eyelidpos = eyelid_traces;

trials.tm = eyelid_times;
trials.fnames = fnames;

trials.c_isi = c_isi;
trials.c_csnum = c_csnum;
trials.c_csdur = c_csdur;
trials.c_csintensity = c_csintensity;
trials.c_usnum = c_usnum;
trials.c_usdur = c_usdur;
trials.c_whitenoise = c_whitenoise;

trials.encoder_displacement = encoder_displacement;
% trials.encoder_time = encoder_time;

trials.trialnum = trialnum;
trials.type = ttype;

trials.iti = diff([0; trialtime]);
trials.iti(trials.trialnum == 1) = 0;

%% baseline correction and normalization %NHE
%baseline correction
traces=bsxfun(@minus, trials.eyelidpos, median(trials.eyelidpos(:,1:40),2));

%normalization and deletion of invalid data
numSTD=3;
tracesnorm=nan(size(traces));

USmean = nanmedian(traces(trials.c_csdur==0, 91:121), 2);
URmean=median(USmean);
tracesnorm = traces/URmean;

stds=nanstd(tracesnorm(:,1:40),[],2); %std of selection above, std of each trial
mus=median(stds); %median value of all the stds
iqrs=iqr(stds); %interquartile range (difference 75th - 25th percentiles of stds).
thresh=mus+numSTD*iqrs; %thresh = median of std + (3*iqrs-difference between interquartile range).
tracesnorm((stds>thresh),:)= nan ;
disp(['Excluded ' num2str(numel(find(stds>thresh))) ' Trials'])
% baseline = mean(mean(trials.eyelidpos(trials.c_csdur==0,1:40)));
% fullclosure = mean(max(trials.eyelidpos(trials.c_csdur==0,96:200),[],2));
% trials.traces = (trials.eyelidpos - baseline) / (fullclosure - baseline);
% trials.meanAll = mean(trials.traces);
% trials.meanCS_US = mean(trials.traces(trials.c_usdur~=0,:));
% trials.meanCSonly = mean(trials.traces(trials.c_usdur==0,:));
% trials.meanUSonly = mean(trials.traces(trials.c_csdur~=0,:));
trials.tracesnorm = tracesnorm;
trials.meanAll = nanmean(trials.tracesnorm);
trials.meanCS_US = nanmean(trials.tracesnorm(trials.c_usdur~=0,:));
trials.meanCSonly = nanmean(trials.tracesnorm(trials.c_usdur==0,:));
trials.meanUSonly = nanmean(trials.tracesnorm(trials.c_csdur==0,:));

%% CR amp calculation
CSselection=tracesnorm(trials.c_csdur==270,75:90);
trials.CRamp=nanmean(tracesnorm(trials.c_csdur==270,75:90),2); %50 ms after CS
meanUR = nanmean(nanmean(tracesnorm(:,101:121),2)) %from 0 to 100 (100 ms aftr US)
trials.CRamp5 = trials.CRamp(trials.CRamp>(0.05*meanUR));

trials.meanCR = nanmean(trials.CRamp)
trials.stdCR = nanstd(trials.CRamp)
trials.meanCR5 = nanmean(trials.CRamp5)
trials.stdCR5 = nanstd(trials.CRamp5)

% ms2frm = @(ms) ms ./ 1e3 * 200;
% pre = 1:ms2frm(200);
% win = ms2frm(300):ms2frm(450);
% trials.CRamp = max(trials.eyelidpos(:, win),[], 2) - mean(trials.eyelidpos(:, pre), 2);
% trials.CRabs = max(trials.eyelidpos(:, win),[], 2);


function db = safestr2double(str)

if ~isempty(str)
    db = str2double(str);
else
    db = 0;
end

