% Reads a timestamp data file and creates a time vector (s).
% (c) Si-yang Yu @ PSI 2018

% Function output: rec_t, defined as an array time points (s)

function rec_t = import_intan_time
  % Get related workspace data from caller
    amplifier_sample_rate = evalin('caller', 'frequency_parameters.amplifier_sample_rate');
  % Execution session
    fileinfo = dir('time.dat');
    num_samples = fileinfo.bytes / 4;
    fid = fopen('time.dat','r');
    rec_t = fread(fid,num_samples,'int32') / amplifier_sample_rate;  % value_type = 'int32'
    fclose(fid);
end