function params = getParams()
% general
params.n_sessions = 2;
params.s = ["S1", "S2"];
params.n_conditions = 2; % There are 2 conditions, CS-US paired and CS-only
params.conditions = [0 1]; % For CS-only, 0, for CS-US paired 1;
params.condition_types = ["CS only", "CS-US combined"];
params.eventTimes = [0 270];
% Paths
params.dirDATA = "/media/mick/DATA/Ilse/spikeSortingUnits/";
params.dirHome = "/home/mick/Desktop/Ilse/spikeSortingUnits/";
params.dataPathDATA = "/media/mick/DATA/Ilse/convertedData/";
params.dataPathHome = "/home/mick/Desktop/Ilse/convertedData/";
params.prbPath = '/home/mick/Desktop/Ilse/DBC_3.1-64-H2_IK.prb';
% Nynke's variables
params.t_pre = 500;
params.t_post = 250; %for CR detection
params.dur = 2800;
params.t_psth = 1000;
params.bin_psth = 150;
params.bin_corr = 250;
% psth and raster variables
params.t_pre_trial = 0.500;
params.t_post_trial = 2;
params.t_raster_edges = 3;
params.raster_spike_height = 0.900;
params.BinW = 0.010;
params.n_trials = 220;
params.t_US = 0.270;
% convert data variables
params.A = 100; %ID number of module % IK: This is the number that is in front of all .continuous files.
params.num_channels = 32;
% mouselist
params.mouseList = ["20-MI19154-08", "20-MI19308-01", "20-MI19601-03", '20-MI19442-03', '20-MI19442-05', '20-MI19442-06', '20-MI19442-08', '21-MI10159-01', '21-MI10159-02', '21-MI10159-06', '21-MI10532-03', '21-MI10532-07', '21-MI16091-04', '21-MI16091-03', '21-MI16183-03', '21-MI16183-05', '22-MI10447-09', '22-MI10447-05', '22-MI10008-08', '22-MI10447-06', '22-MI11756-06', '22-MI11756-07', '22-MI12417-08', '22-MI12410-07', '22-MI13134-01', '22-MI13134-03', '22-MI13989-05', '22-MI13989-07', '22-MI14020-03', '22-MI14020-07', '22-MI14020-08', 'MI22.02438.04', 'MI22.02438.05', 'MI23.00102.01', 'MI23.00102.03', 'MI22.01127.02', 'MI22.01127.03', 'MI22.01127.04', 'MI22.01125.07', 'MI22.02178.02', 'MI22.02178.01', 'MI22.02178.04', 'MI22.02178.03', 'MI22.02170.01', 'MI23.00244.04', 'MI23.01047.03', 'MI23.01712.07', 'MI23.01753.03', 'MI23.01412.04', 'MI23.01753.02', 'MI23.02063.04', 'MI23.02063.05'];
params.mouseNames = ["Aurora", 'Bacchus', 'Cupid', 'Diana.1', 'Diana.2', 'Diana.3', 'Diana.4', 'Epona.1', 'Epona.2', 'Epona.3', 'Fortuna.1', 'Fortuna.2', 'Genius.2', 'Genius.1', 'Hercules.1', 'Hercules.2', 'Invidia.3', 'Invidia.1', 'Juno', 'Invidia.2', 'Luna.1', 'Luna.2', 'Mars', 'Neptune', 'Orcus.1', 'Orcus.2', 'Pluto.1', 'Pluto.2', 'Quiritis.1', 'Quiritis.2', 'Quiritis.3', 'Roma.1', 'Roma.2', 'Saturn.1', 'Saturn.2', 'Trivia.1', 'Trivia.2', 'Trivia.3', 'Trivia.4', 'Vulcan.2', 'Vulcan.1', 'Vulcan.4', 'Vulcan.3', 'Ares', 'Bia', 'Chronos', 'Demeter', 'Eros.2', 'Fury', 'Eros.1', 'Gaia.1', 'Gaia.2'];
end
