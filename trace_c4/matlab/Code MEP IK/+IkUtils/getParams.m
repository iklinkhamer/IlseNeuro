%% IK 8-4-24
function params = getParams(kwargs)
arguments
    % General settings for K-means clustering and filtering
    kwargs.Sspk_base_filter (1,1) double = 20;
    kwargs.Sspk_max_filter (1,1) double = 200; % IK note: this might be too low. Check
    kwargs.validTrialThres = 80; % Trials that are at least 20Hz are valid.
    
    % General settings for histogram plots
    kwargs.histogramModulationThreshold_sspks = 5; % 5 Hz taken from ten Brinke et al., 2017 % IK added. Check if this value is correct later.
    kwargs.minimumModulationTimeThreshold_sspks = 10; % 10 ms taken from ten Brinke et al., 2017 % IK added.
    kwargs.sspkRanges = struct ...
        ( cs = struct ...
            ( min = 0.020 ... % IK changed from 0.055 and from 0.05
            , max = 0.250 ...
            ) ...
        , us = struct ...
            ( min = 0.250... % IK change %0.05 ...
            , max = 0.350... %0.150... % IK change. %0.070 ... % I.e. 50 milliseconds after the US event
            ) ...
        );
    
    
    % Trace event times
    kwargs.delayEventsCs = [0, 0.25];  % IK change
    kwargs.delayEventsUs = [-0.25, 0]; % IK change
    
    
    % Settings for histogram bin width
    kwargs.BinW = 0.005; % 5 ms
    kwargs.sd_thres = 3/4 * std(rand(1e6,1)); % std for uniformly distributed spiking events % IK note: check for simple
        
    % Params for contamination
    kwargs.cThreshold = 10;
    kwargs.cBinW = 5;
    kwargs.cBinLimits = [0,1200];
    
end

    kwargs.psthRanges = struct ...
        ( cs_full = struct ...
            ( min = -0.25 ...
            , max = 0.75... %0.75 ... % IK change
            ) ...
        , cs_zoomed = ...
            struct ...
                ( min = min(-0.02   , kwargs.sspkRanges.cs.min) ...
                , max = max(0.5 , kwargs.sspkRanges.cs.max) ... % IK change
                ) ...
        , us_full = ...
            struct ...
                ( min = -0.75 ...
                , max = 0.25 ...
                ) ...
        , us_zoomed = ...
            struct ...
                ( min = -.02 ...
                , max = 0.200 ...
                ) ...
        );
%% overviewTraces
kwargs.lastDay = 10;
kwargs.nTimeSteps = 200;
kwargs.nTrials = 240;
kwargs.nCSplusRegTrials = 220;

kwargs.CSdur = 270;
kwargs.USdur = 20;

kwargs.thresCRperc = 5;
kwargs.thresCRonset = 5;

kwargs.tracesRanges = struct ...
        ( baseline = struct ...
            ( min = 1 ... 
            , max = 40 ...
            ) ...
        , cs = struct ...
            ( min = 41 ... 
            , max = 90 ...
            ) ...
        , us = struct ...
            ( min = 91 ... 
            , max = 140 ...  
            ) ...
         , after = struct ...
            ( min = 141 ...
            , max = 200 ...
            ) ...
        );

% kwargs.complex_edges = [0:kwargs.BinW:0.2]; % this no longer works since we
% have multiple complex spike search ranges (cs, prior, us).

% general
kwargs.n_sessions = 3;
kwargs.s = ["S1", "S2", "S3"];
kwargs.n_conditions = 3; % There are 3 conditions, CS-US paired, CS-only and US-only
kwargs.conditions = [0 1 2]; % For CS-only, 0, for CS-US paired 1, for US-only 2;
kwargs.condition_types = ["CS only", "CS-US combined", "US only"];
kwargs.eventTimes = [0 250];
% Paths
kwargs.dirDATA = "/mnt/Data/Ilse/";
kwargs.dirHome = "/home/i.klinkhamer/Documents/";
kwargs.dirHomeData = fullfile(kwargs.dirHome, "Data");
kwargs.dirDATAData = fullfile(kwargs.dirDATA, "Data");
kwargs.dirHomeCode = fullfile(kwargs.dirHome, "ephys_code");
kwargs.pathSpikeSortingDATA = fullfile(kwargs.dirDATAData, "spikeSortingUnits");
kwargs.pathSpikeSortingHome = fullfile(kwargs.dirHomeData, "spikeSortingUnits");
kwargs.pathEphysDataDATA = fullfile(kwargs.dirDATAData, "convertedData");
kwargs.pathEphysDataHome = fullfile(kwargs.dirHomeData, "ephysData");
kwargs.pathBehaviorDataDATA = fullfile(kwargs.dirDATAData, "behaviorData");
kwargs.pathBehaviorDataHome = fullfile(kwargs.dirHomeData, "behaviorDataResults");
kwargs.figPath = fullfile(kwargs.dirHome, "Figures");
kwargs.prbPath = fullfile(kwargs.dirHome, "DBC_3.1-64-H2_IK.prb");
% Nynke's variables
kwargs.t_pre = 500;
kwargs.t_post = 250; %for CR detection
kwargs.dur = 2800;
kwargs.t_psth = 1000;
kwargs.bin_psth = 150;
kwargs.bin_corr = 250;
kwargs.isi = 250;
kwargs.digital_US_delay_ms = 12; % US appears 12 ms earlier in the digitally recorded US times than it appears in real life.
kwargs.digital_US_delay = kwargs.digital_US_delay_ms*0.001;
kwargs.CS = 1;
kwargs.US = 3;
% psth and raster variables
kwargs.t_pre_trial = 0.500;
kwargs.t_post_trial = 2;
kwargs.t_raster_edges = 3;
kwargs.t_US_offset = 0.250;
kwargs.raster_spike_height = 0.900;
kwargs.histogramHeight = 15;
% kwargs.BinW = 0.010;
kwargs.t_US = 0.250;
kwargs.spikeTimesBinW = 0.0005;
% convert data variables
kwargs.num_channels = 32;
% mouselist
kwargs.mouseList = ["20-MI19154-08", "20-MI19308-01", "20-MI19601-03", '20-MI19442-03', '20-MI19442-05', '20-MI19442-06', '20-MI19442-08', '21-MI10159-01', '21-MI10159-02', '21-MI10159-06', '21-MI10532-03', '21-MI10532-07', '21-MI16091-04', '21-MI16091-03', '21-MI16183-03', '21-MI16183-05', '22-MI10447-09', '22-MI10447-05', '22-MI10008-08', '22-MI10447-06', '22-MI11756-06', '22-MI11756-07', '22-MI12417-08', '22-MI12410-07', '22-MI13134-01', '22-MI13134-03', '22-MI13989-04','22-MI13989-05', '22-MI13989-07', '22-MI14020-03', '22-MI14020-07', '22-MI14020-08', 'MI22.02438.04', 'MI22.02438.05', 'MI23.00102.01', 'MI23.00102.03', 'MI22.01127.02', 'MI22.01127.03', 'MI22.01127.04', 'MI22.01125.07', 'MI22.02178.02', 'MI22.02178.01', 'MI22.02178.04', 'MI22.02178.03', 'MI22.02170.01', 'MI23.00244.04', 'MI23.01047.03', 'MI23.01712.07', 'MI23.01753.03', 'MI23.01412.04', 'MI23.01753.02', 'MI23.02063.04', 'MI23.02063.05', 'MI23.02436.01', 'MI23.02436.03', 'MI23.02436.04', 'MI23.03123.02', 'MI23.03123.03', 'MI23.03123.04', 'MI23.03374.01', 'MI23.03374.03', 'MI23.05990.02', 'MI23.05990.03', 'MI23.05990.04', 'MI23.05990.05', 'MI23.05990.08'];
kwargs.mouseNames = ["Aurora", 'Bacchus', 'Cupid', 'Diana.1', 'Diana.2', 'Diana.3', 'Diana.4', 'Epona.1', 'Epona.2', 'Epona.3', 'Fortuna.1', 'Fortuna.2', 'Genius.2', 'Genius.1', 'Hercules.1', 'Hercules.2', 'Invidia.3', 'Invidia.1', 'Juno', 'Invidia.2', 'Luna.1', 'Luna.2', 'Mars', 'Neptune', 'Orcus.1', 'Orcus.2', 'Pluto.1', 'Pluto.2', 'Pluto.3', 'Quiritis.1', 'Quiritis.2', 'Quiritis.3', 'Roma.1', 'Roma.2', 'Saturn.1', 'Saturn.2', 'Trivia.1', 'Trivia.2', 'Trivia.3', 'Trivia.4', 'Vulcan.2', 'Vulcan.1', 'Vulcan.4', 'Vulcan.3', 'Ares', 'Bia', 'Chronos', 'Demeter', 'Eros.2', 'Fury', 'Eros.1', 'Gaia.1', 'Gaia.2', 'Hades.1', 'Hades.2', 'Hades.3', 'Icarus.1', 'Icarus.2', 'Icarus.3', 'Jason.1', 'Jason.2', 'Kratos.1', 'Kratos.2', 'Kratos.3', 'Kratos.4', 'Kratos.5'];
kwargs.Shank2MUT = ["Diana.1", 'Diana.2', 'Diana.4', 'Epona.1', 'Epona.2', 'Fortuna.2', 'Genius.2', 'Genius.1', 'Hercules.2', 'Invidia.3', 'Juno', 'Luna.1', 'Orcus.1', 'Pluto.2', 'Quiritis.2', 'Quiritis.3', 'Saturn.2', 'Bia', 'Fury'];
kwargs.Shank2WT = ["Diana.3", 'Epona.3', 'Fortuna.1', 'Hercules.1', 'Invidia.1', 'Invidia.2', 'Luna.2', 'Orcus.2', 'Pluto.1', 'Quiritis.1', 'Saturn.1', 'Chronos', 'Demeter'];

types = struct;
types.full = ["Shank2KOWT", "Shank2KOMUT", "Shank2L7WT", "Shank2L7MUT", "Tsc1KOWT", "Tsc1KOMUT", "Tsc1L7WT", "Tsc1L7MUT"];
types.gene = ["Shank2", "Shank2", "Shank2", "Shank2", "Tsc1", "Tsc1", "Tsc1", "Tsc1"];
types.loc = ["KO", "KO", "L7", "L7", "KO", "KO", "L7", "L7"];
types.type = ["WT", "MUT", "WT", "MUT", "WT", "MUT", "WT", "MUT"];
kwargs.types = types;

kwargs.openEphysMapping = [50, 49, 53, 55, 57, 59, 62, 64, 63, 61, 60, 58, 56, 54, 51, 52, 16, 15, 11, 9, 7, 5, 4, 2, 1, 3, 6, 8, 10, 12, 13, 14];

params = kwargs;

end