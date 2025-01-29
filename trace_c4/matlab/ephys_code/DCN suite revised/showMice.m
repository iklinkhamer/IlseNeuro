% Show registered mice and status of their data curations.

% Copyright 2023 NarainLab
% Erasmus Medical Center, Rotterdam, The Netherlands
% v2023-01-04
function showMice()
    
    mouseList = defaultMice();
    
    nMice = numel(mouseList);
    
    IkUtils.printIndented ...
        ( sprintf("Found %d registered mice:\n", nMice) ...
        )
    
    arrayfun(@showMouse, mouseList);
    
end