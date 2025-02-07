function flattened = flatten(cell_array)
    flattened = cell(0);
    for i = 1:numel(cell_array)
        if iscell(cell_array{i})
            flattened = [flattened, flatten(cell_array{i})];
        else
            flattened = [flattened, cell_array{i}];
        end
    end
end
