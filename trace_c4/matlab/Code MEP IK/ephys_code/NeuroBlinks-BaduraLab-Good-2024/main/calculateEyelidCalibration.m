function calculateEyelidCalibration()
% -- load data again --
metadata = getappdata(0, 'calibration_metadata');
vid = getappdata(0, 'calibration_vid');
vid_ts = getappdata(0, 'calibration_vid_ts');

% --- eyelid trace --
[trace, t] = vid2eyetrace(vid, metadata, metadata.cam(1).thresh);

% --- cal ---
ind_t = t < 0.1;
calib_offset = min(trace(ind_t));
maxclosure = max(trace(ind_t));
calib_scale = maxclosure - calib_offset;

% -- save cal data to root metadata ---
metadata.cam(1).calib_offset = calib_offset;  
metadata.cam(1).calib_scale = calib_scale;
metadata.cam(1).cal = 0;    % Clear calibration flag

setappdata(0,'metadata',metadata);

fprintf('calib_offset = %d.  calib_scale = %d.\n', calib_offset, calib_scale)
fprintf('thresh = %d.\n', round(metadata.cam(1).thresh * 255))

videoname = sprintf('%s\\%s_calib', metadata.folder, metadata.basename);
save(videoname, 'vid', 'vid_ts', 'metadata', '-v6');    

fprintf('Video from calibration trial successfully written to disk.\n')