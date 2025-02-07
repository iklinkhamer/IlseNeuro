% Import JRCLUST generated CSV file for further analysis.
% (c) Si-yang Yu @ PSI 2018

% Function input:  csv_file, defined as *.csv file generated by JRCLUST
% Function output: cell_spk, defined as an array of struct (3 fields, defined below)
%                              nr: the cell number, 0 means noise channel
%                              t:  the spike timings
%                              ch: the channel(s) where the cell was recorded

function cell_spk = import_jrc_csv(csv_file)
  % Read in file and prepare data set
    spk_all = csvread(csv_file);
    n_cell = max(spk_all(:,2));
  % Related data extraction loop
    i = 0;
    cell_spk(n_cell + 1) = struct('nr',[],'t',[],'ch',[]);
    while i <= n_cell
        t_idx = spk_all(:,2) == i;
        cell_t = spk_all(t_idx,1);
        cell_ch = unique(spk_all(t_idx,3));
        cell_spk(i+1) = struct('nr',i,'t',cell_t,'ch',cell_ch);
        i = i + 1;
    end
end