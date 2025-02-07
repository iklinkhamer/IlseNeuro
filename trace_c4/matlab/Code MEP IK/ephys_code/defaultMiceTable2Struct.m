table_path = "/home/i.klinkhamer/Documents/Data/defaultMiceStruct.csv";
opts = detectImportOptions(table_path);

% Treat all columns as character vectors
opts.VariableTypes = repmat({'char'}, 1, numel(opts.VariableNames));

T = readtable(table_path, opts);
% Assuming T is your table read from the file

% Initialize an empty cell array to store the ephys dates
ephysdates = cell(height(T), 1);

% Iterate over each variable in the table
varNames = T.Properties.VariableNames;
for i = 1:numel(varNames)
    varName = varNames{i};
    % Check if the variable name contains 'ephysdate'
    if contains(varName, 'ephysdate', 'IgnoreCase', true)
        % Concatenate the data into the ephysdates cell array
        ephysdates(:, end+1) = table2cell(T(:, varName));
        T = removevars(T, varName);
    end
end

% Remove the first empty column
ephysdates = ephysdates(:, 2:end);

% Add ephysdates as a new field to the table
T.ephysdates = ephysdates;

S = table2struct(T);

 allEmpty = all(cellfun(@isempty, S(1).ephysdates));