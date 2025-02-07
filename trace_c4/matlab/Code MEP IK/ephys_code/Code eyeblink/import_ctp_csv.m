% Import cell type CSV file for further analysis.
% (c) Si-yang Yu @ PSI 2018

% Function input:  csv_file, defined as *.csv file containing cell type information
% Function output: cell_tp,  defined as an array of struct (2 fields, defined below)
%                              nr:   the cell number
%                              type: the cell type

function cell_tp = import_ctp_csv(csv_file)
    table_in = readtable(csv_file,'ReadVariableNames',0);
    table_in.Properties.VariableNames = {'nr','type'};
    cell_tp = transpose(table2struct(table_in));
end