clear all
close all

%% load data
directory='C:\Users\Badura blink\Documents\Open Ephys\test1_2020-11-10_16-04-05\Record Node 106\experiment1\recording1\spikes\Spike_Sorter-104.0\spike_group_1';
metadata=int64(readNPY(strcat(directory,'\metadata.npy')));
spike_clusters=int64(readNPY(strcat(directory,'\spike_clusters.npy')));
spike_electrode_indices=int64(readNPY(strcat(directory,'\spike_electrode_indices.npy')));
spike_times=int64(readNPY(strcat(directory,'\spike_times.npy')));
spike_waveforms=int64(readNPY(strcat(directory,'\spike_waveforms.npy')));

data=array2table([metadata spike_clusters spike_electrode_indices spike_times spike_waveforms(:,:,1)]);
data.Properties.VariableNames{1} = 'metadata';
data.Properties.VariableNames{4} = 'clusters';
data.Properties.VariableNames{5} = 'indices';
data.Properties.VariableNames{6} = 'time';
data.Properties.VariableNames{7} = 'waveform1';
data.Properties.VariableNames{8} = 'waveform2';

%% plot
ind=linspace(1,15,8);
for i=1:8   
    for j=1:2
        subplot(4,4,(ind(i)-1+j))
        if j==1
            plot(data.time(data.indices==i,:), data.waveform1(data.indices==i,:));
        else
            plot(data.time(data.indices==i,:), data.waveform2(data.indices==i,:));
        end
        xlim([1.5e7 1.505e7]);
    end
end