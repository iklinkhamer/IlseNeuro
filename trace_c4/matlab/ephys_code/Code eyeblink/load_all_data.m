%% load all data
clear all 
close all

num_channels=32;
A=100; %ID number of module
%B=2; %experiment number, only necessary for multiple runs with same name
path='D:\Nynke\Ephys-data\MI19442-06F_2021-05-31_13-28-33_2\Record Node 101\'; %insert path to files here
formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
%formatSpec=''\???.spike'; %filename format for .spike data
chSel = [4 24 27];

for i=1:3
    ch=chSel(i)
    filepath=sprintf(formatSpec,path,A,(ch)) %generate file path
    [data(i,:), timestamps(i,:), info(i,:)]=load_open_ephys_data_chunked(filepath, 0, 1800); %load data
end

%% plotjesss
figure ()
title('Signal vs time(s)')
for i=1:32
    subplot(4,round(num_channels/4),i)
    plot(timestamps(i,1:30000),data(i,1:30000), 'color', rand(1,3))
    ylim([-1000 1000])
    title(i)
    %xlabel('time(s)')
end


%% saving
%save(fullfile('D:\Nynke\Ephys-data', 'MI19442-03F_ephys1.dat'), 'data', '-v7.3');
% fid = fopen('YourOutput.dat', 'w');
% fwrite(fid, data);
% fclose(fid);
% 
%  save mymatrix.dat data -ascii

fileID = fopen(strcat(path, 'RAWdata.bin'),'w');
fwrite(fileID,data, 'double');
fclose(fileID);

% %%
% fileID = fopen('datatest.bin');
% datatest = fread(fileID, 'double');

%%

[data2, timestamps2, info2]=load_open_ephys_data(strcat(path, 'all_channels.events'));

times = timestamps2 - timestamps(1,1);

%%
save(strcat(path,'stimtimes.mat'),'times')

T=timestamps(:,1);
save(strcat(path, 'timestamps.mat'), 'T')