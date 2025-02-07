%% load all data
clear all

num_channels=32;
A=100; %ID number of module
%B=2; %experiment number, only necessary for multiple runs with same name
 %insert path to files here
prompt = "Which mouse's data do you want to convert? ";
mname = input(prompt,"s");
prompt2 = "Which session do you want to convert? ";
session = input(prompt2);

path = getMousePath(mname, session);
%path =['/home/mick/Desktop/10008-08F_2022-05-23_11-05-23_1/Record/'];
% try
%     formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
% catch
%     formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
% end

chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];

%chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];


% try
    formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
    for i=1:32
        ch=chSel(i)
        filepath=sprintf(formatSpec,path,A,(ch)) %generate file path
        [data(i,:), timestamps(i,:), info(i,:)]=load_open_ephys_data(filepath); %load data
    end
% catch
%     formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
%     for i=1:32
%         ch=chSel(i)
%         filepath=sprintf(formatSpec,path,A,(ch)) %generate file path
%         [data(i,:), timestamps(i,:), info(i,:)]=load_open_ephys_data(filepath); %load data
%     end
% end

%% preview data
figure ()
title('Signal vs time(s)')
for i=1:32
    subplot(8,round(num_channels/8),i)
    plot(timestamps(i,1:15000),data(i,1:15000),'color', 'black')
    ylim([-400 400])
    title(i)
    %xlabel('time(s)')
end



fileID = fopen(strcat(path, 'RAWdata.bin'),'w');
fwrite(fileID,data, ['double']);
fclose(fileID);

% %%
% fileID = fopen('datatest.bin');
% datatest = fread(fileID, 'double');

%%

[data2,timestamps2, info2]=load_open_ephys_data(strcat(path, 'all_channels.events'));

times = timestamps2 - timestamps(1,1);

%%
save(strcat(path,'stimtimes.mat'),'times')

T=timestamps(:,1);
save(strcat(path, 'timestamps.mat'), 'T')
