clear all
close all

%% continuous 
directoryC='C:\Users\Badura blink\Documents\Open Ephys\MI134251-08_2020-12-18_11-17-49\Record Node 101\experiment1\recording1\continuous\Rhythm_FPGA-100.0';
[dataCont, timestamps, syncTimestamps]=readContinuous(directoryC);
plotContinuous(timestamps, dataCont)

%% continuous after BPfilter
directoryC='C:\Users\Badura blink\Documents\Open Ephys\MI134251-08_2020-12-18_11-17-49\Record Node 106\experiment1\recording1\continuous\Rhythm_FPGA-100.0';
[dataContBP, timestamps, syncTimestamps]=readContinuous(directoryC);
plotContinuous(timestamps, dataContBP)

% %% spikes
% directory='I:\Nynke\Open Ephys\test5_2020-11-26_10-43-08\Record Node 106\experiment1\recording1\spikes\Spike_Sorter-104.0\spike_group_1';
% metadata=int64(readNPY(strcat(directory,'\metadata.npy')));
% spike_clusters=int64(readNPY(strcat(directory,'\spike_clusters.npy')));
% spike_electrode_indices=int64(readNPY(strcat(directory,'\spike_electrode_indices.npy')));
% spike_times=int64(readNPY(strcat(directory,'\spike_times.npy')));
% spike_waveforms=int64(readNPY(strcat(directory,'\spike_waveforms.npy')));
% 
% dataS=array2table([metadata spike_clusters spike_electrode_indices spike_times spike_waveforms(:,:,1)]);
% dataS.Properties.VariableNames{1} = 'metadata';
% dataS.Properties.VariableNames{4} = 'clusters';
% dataS.Properties.VariableNames{5} = 'indices';
% dataS.Properties.VariableNames{6} = 'time';
% dataS.Properties.VariableNames{7} = 'waveform';

% plotSpikes(dataS)

%% Functions
function [dataCont, timestamps, syncTimestamps]=readContinuous(directory)
timestamps=readNPY(strcat(directory, '\timestamps.npy'));
syncTimestamps=readNPY(strcat(directory, '\synchronized_timestamps.npy'));
fileID=fopen(strcat(directory, '\continuous.dat'));
dataCont=fread(fileID, [64, length(timestamps)],'int64');
end

function plotContinuous(timestamps, dataCont)
timeCont=timestamps;%double(timestamps)/30e3;
figure
for i=1:64 %size(dataCont2,1)
subplot(8,8,i)
plot(1:length(dataCont), dataCont(i,:))
ylim([-max(max(dataCont))-100 max(max(dataCont))+100])
ylabel(i)
xlim([10000 20000])
end

figure
plot(timeCont(1:length(dataCont)), dataCont(1:20,:))

end

function plotSpikes(dataS)
figure
for i=1:64
    subplot(8,8,i)
    plot(dataS.time(dataS.indices==i,:), dataS.waveform(dataS.indices==i,:));
    ylabel(i)
end
end
