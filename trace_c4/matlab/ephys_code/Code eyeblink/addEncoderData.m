function trials = addEncoderData(trials, varargin)
% trials = addEncoderData(trials, {directory})
% Add encoder data to trials struct using either same directory as fnames in trial struct or optionally supplied directory to locate encoder MAT files

% Check if we need to use a different base directory
if length(varargin) > 0
    use_original_path = 0;
    base_path = varargin{1};
else
    use_original_path = 1;
end

fnames = trials.fnames;

trials.encoder_displacement = NaN(size(trials.eyelidpos));

for i=1:length(fnames)

    [p,fname,ext] = fileparts(fnames{i});

    if use_original_path
        base_path = p;
    end

    basename = regexp(fname, '[A-Z]\d\d\d_\d\d\d\d\d\d_[ts]\d\d[a-z]?_\d\d\d', 'match', 'once');

    encoder_fname = fullfile(base_path, sprintf('%s_encoder.mat', basename));

    if exist(encoder_fname, 'file')

        load(encoder_fname);

        trials.encoder_displacement(i, :) = encoder.displacement;

    end

end