%% Mick & Ilse 7/9/23     load all data
function convertData(kwargs)
arguments
%     kwargs.session = 1;
%     kwargs.firstTimeSorting = true;
    kwargs.oldData = true;
    kwargs.num_channels=32;
    kwargs.A=100; %ID number of module
end

if kwargs.oldData == 1
    prompt = "Enter mouse index: ";
    mname = input(prompt, 's');
    prompt2 = "Enter session number: ";
    session = input(prompt2);

    directory = "/home/mick/Desktop/spikeSortedUnits/";

    if session == 1
        session_n = "01";
    elseif session == 2
        session_n = "02";
    end

    if isfile(directory + mname + "/" + session_n + "/" + "100_CH1.continuous")
        formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
        chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
    elseif isfile("Desktop/spikeSortedUnits/" + mname + "/" + session_n + "/" + "continuous.dat")
        formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
        chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
    end
    %
%     for i = 1:length(delimiters)
%         if dataType == "dat"
%             filetest = directory + name_idx_1 + delimiters(i) + name_idx_2 + "/" + session_n + "/" + "Config_h32_oe_dat.prm";
%         elseif dataType == "bin"
%             filetest = directory + name_idx_1 + delimiters(i) + name_idx_2 + "/" + session_n + "/" + "Config_h32_oe.prm";
%         end
%         %f = dir(filetest);
%         if isfile(filetest)
%             filename = filetest;
%         end
%     end

    if isfile(erase(filename, ".prm") + "_full.prm")
        str = "irc manual " + filename;
    elseif kwargs.firstTimeSorting == 0
        str = "irc full " + filename;
    end
    



%B=2; %experiment number, only necessary for multiple runs with same name
 %insert path to files here
prompt = "Which mouse's data do you want to convert? ";
mname = input(prompt,"s");
prompt2 = "Which session do you want to convert? ";
session = input(prompt2);


[path,chSel_type, formatSpec_type] = getMousePath(mname, session, whose);

if formatSpec_type == 1
     formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
elseif formatSpec_type ==2
     formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
end

if chSel_type == 32
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
elseif chSel_type == 64
    chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
end


for i=1:16
    ch=chSel(i);
    filepath=sprintf(formatSpec,path,kwargs.A,(ch)); %generate file path
    [data1(i,:), timestamps1(i,:), ~]=load_open_ephys_data(filepath); %load data
end
save(strcat(path, 'data1.mat'),'data1', '-v7.3'); 
save(strcat(path, 'timestamps1.mat'),'timestamps1', '-v7.3'); 
data1 = [];
timestamps1 = [];
for i=1:16
    ii = i + 16;
    ch=chSel(ii);
    filepath=sprintf(formatSpec,path,A,(ch)); %generate file path
    [data2(i,:), timestamps2(i,:), ~]=load_open_ephys_data(filepath); %load data
end
save(strcat(path, 'data2.mat'),'data2', '-v7.3'); 
save(strcat(path, 'timestamps2.mat'),'timestamps2', '-v7.3'); 
timestamps2 = [];

%%
%path ='/home/mick/Desktop/10008-08F_2022-05-23_11-05-23_1/Record/';
data1 = load(strcat(path, 'data1.mat'));
data1 = data1.data1;
% data2 = load(strcat(path, 'data2.mat'));
% data2 = data2.data2;
data = [data1;data2];
save(strcat(path, 'data.mat'),'data', '-v7.3'); 
data1 = [];
data2 = [];
fileID = fopen(strcat(path, 'RAWdata.bin'),'w');
fwrite(fileID,data, 'double');
fclose(fileID);

timestamps1 = load(strcat(path, 'timestamps1'));
timestamps1 = timestamps1.timestamps1;
timestamps2 = load(strcat(path, 'timestamps2'));
timestamps2 = timestamps2.timestamps2;
timestamps = [timestamps1;timestamps2];
timestamps1 = [];
timestamps2 = [];
save(strcat(path, 'timestampstest.mat'),'timestamps', '-v7.3'); 
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



% %%
% fileID = fopen('datatest.bin');
% datatest = fread(fileID, 'double');

%%
path ='/home/mick/Desktop/10008-08F_2022-05-23_11-05-23_1/Record/';
  chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
  A = 100;
  formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
ch=chSel(1);
filepath=sprintf(formatSpec,path,A,(ch)); %generate file path
[~,timestamps(1,:), ~]=load_open_ephys_data(filepath); %load data

[data2,timestamps2, info2]=load_open_ephys_data(strcat(path, 'all_channels.events'));

times = timestamps2 - timestamps(1,1);

%%
save(strcat(path,'stimtimes.mat'),'times')

T=timestamps(:,1);
save(strcat(path, 'timestamps.mat'), 'T')
