%% Ilse 8/9/23   convert events to .mat files
function openephys_events_only_IK(kwargs)
arguments
    kwargs.whose = "Ilse";
    kwargs.num_channels = 32;
    kwargs.A = 100; %ID number of module
    kwargs.previewData = false;
    kwargs.saveConvertedData = true;
    kwargs.saveStimTimestamps = true;
    kwargs.oldData = true;
end


if kwargs.oldData == 1
    prompt = "Enter mouse index: ";
    mname = input(prompt, 's');
    prompt2 = "Enter session number: ";
    session = input(prompt2);

    directory = '/home/mick/Desktop/spikeSortedUnits/';
   
    if session == 1
        session_n = 'Session_1';
    elseif session == 2
        session_n = 'Session_2';
    end

    path = strcat(directory,mname,'/',session_n,'/'); 

    if isfile(path + "100_CH1.continuous")
        formatSpec='%s%d_CH%d.continuous'; %filename format for .continuous data
        chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
    elseif isfile(path + "100_1.continuous")
        formatSpec='%s%d_%d.continuous'; %filename format for .continuous data
        chSel = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32];
    end
    
[~,timestamps2, ~]=load_open_ephys_data(strcat(path, 'all_channels.events'));

filepath=sprintf(formatSpec,path,kwargs.A,(1)); %generate file path
[~,timestamps(1,:), ~]=load_open_ephys_data(filepath); %load data

times = timestamps2 - timestamps(1,1);

save(strcat(path,'stimtimes.mat'),'times')

T=timestamps(:,1);
save(strcat(path, 'timestamps.mat'), 'T')
end