clear;

% Input parameters
%   csv_name = 'testrun.csv';
  t_pre = 500;
  t_post = 250; %for CR detection
  dur = 2000;
  t_psth = 1000;
  bin_psth = 150;
  bin_corr = 250;

% Import Intan header and path
  read_Intan_RHD2000_file;
  cd(path);
% Import Intan time vectors
  t = import_intan_time;
% Import Intan ADC channels vector
  adc = import_intan_adc([]);
% Import CSV file from JRCLUST
  spk = import_jrc_csv(csv_name);
  ctp = import_ctp_csv('cell_tp.csv');
% Detection CS/US edge
  cs_lc = find_trg_pks([adc(1).v],3,t,dur);
  us_lc = find_trg_pks([adc(2).v],3,t,dur);
  us_lc(10:10:end,:) = [];
% Detect CR/UR in eyeblink trace
  blk = blk_detn(v,t,cs_lc,us_lc,t_pre,t_post,350,dur); %-500 to 2000ms trace
  %blk_plot(t,blk);
% Calculte non_CR trial PSTH (not including CR_only trials), CR trial PSTH (not including CR_only trials) and CS_only trial PSTH
  [non_cr_locs,cs_only_locs,Clocs] = locs(blk,cs_lc);
  [Cpsth,Ctas] = psth_calc(spk,Clocs,50,5,t_pre,t_psth);
  [NCpsth,NCtas] = psth_calc(spk,non_cr_locs,50,5,t_pre,t_psth);
  [CSpsth,CStas] = psth_calc(spk,cs_only_locs,50,5,t_pre,t_psth);
