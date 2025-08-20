filename = fullfile(Env.getBayesLabUserRoot, 'ContextMouseExperiments', 'EphysSessions', 'L_test_2025-07-31_13-43-30_2', 'Record Node 121', 'experiment1', 'recording1', 'continuous', 'Acquisition_Board-108.acquisition_board', 'continuous.dat');
num_channels = 32;
sampling_rate = 30000;
close all
fid = fopen(filename, 'r');
raw_data = fread(fid, [num_channels, Inf], 'int16');
fclose(fid);

raw_data = raw_data';  % samples x channels

start_sample = 5;
duration_seconds = 2;
end_sample = start_sample + duration_seconds * sampling_rate - 1;

time_vector = (start_sample:end_sample) / sampling_rate;

channels_to_plot = 1:32;
offset_step = 50000;

figure('Position', [100, 100, 1600, 400]);

% Raw data plot (no filter)
subplot(2,1,1)
hold on
for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    raw_ch = double(raw_data(start_sample:end_sample, ch));
    plot(time_vector, raw_ch + (i-1)*offset_step)
end
xlabel('Time (s)')
ylabel('Amplitude + offset')
title('Raw data, no normalization or filtering')
yticks([])
hold off

% Filtered data plot (300-6000 Hz)
subplot(2,1,2)
[b,a] = butter(3, [300, 6000]/(sampling_rate/2), 'bandpass');
hold on
for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    raw_ch = double(raw_data(start_sample:end_sample, ch));
    filt_ch = filtfilt(b, a, raw_ch);
    filt_ch = (filt_ch - mean(filt_ch)) / std(filt_ch);
    plot(time_vector, filt_ch + (i-1)*offset_step)
end
xlabel('Time (s)')
ylabel('Amplitude + offset')
title('Filtered data (300–6000 Hz)')
yticks([])
hold off

% PSD of one channel
ch = 1;
raw_ch = double(raw_data(start_sample:end_sample, ch));
figure;
pwelch(raw_ch, [], [], [], sampling_rate);
title('PSD of raw channel data');

% New figure: highlight parts of filtered signal above noise
% Parameters for spike window
pre_samples = round(0.5e-3 * sampling_rate);   % 0.5 ms before threshold crossing
post_samples = round(0.5e-3 * sampling_rate);  % 0.5 ms after

figure('Position', [100, 600, 1600, 400]);
hold on

for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    raw_ch = double(raw_data(start_sample:end_sample, ch));
    filt_ch = filtfilt(b, a, raw_ch);
    filt_ch = filt_ch - mean(filt_ch);

    noise_std = std(filt_ch);
    threshold = 3 * noise_std;

    % Find threshold crossings (positive and negative)
    crossings = find(abs(filt_ch) > threshold);

    % Skip close ones (prevent overlap)
    refractory = round(1e-3 * sampling_rate);  % 1 ms refractory
    if numel(crossings) > 1
        crossings = crossings([true; diff(crossings) > refractory]);
    end


    % Plot all spikes
    for j = 1:length(crossings)
        idx = crossings(j);
        if idx - pre_samples < 1 || idx + post_samples > length(filt_ch)
            continue
        end
        window = idx - pre_samples : idx + post_samples;
        t_win = time_vector(window);
        spike = filt_ch(window) + (i-1)*offset_step;
        plot(t_win, spike, 'r')
    end
end

xlabel('Time (s)')
ylabel('Amplitude + offset')
title('Full spike waveforms above 3× noise level (300–6000 Hz)')
yticks([])
hold off

