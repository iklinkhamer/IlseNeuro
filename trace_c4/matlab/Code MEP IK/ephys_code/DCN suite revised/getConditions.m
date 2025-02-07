% Julius: this function does more than its name specifies
% TODO: rename this function or splice out the different functionalities into
% separate functions.
function [nConditions, index, distType, histColor] = getConditions(unit)

try
    nConditions = size(unit.neuron(1).RasterXY_cs,1);
catch
    nConditions = 2;
end
% Delta mouse or Uniform mouse

% switch nConditions
%     case 2
%         distType = "Delta";
%         histColor = getColors().c_delta;
%     case 6
%         distType = "Uniform";
%         histColor = getColors().c_uniform;
%     otherwise
%         error ...
%             ( "Could not assign distribution type: unknown number of conditions: %d" ...
%             , nConditions ...
%             )
% end

% if (nConditions > 2), distType = "Uniform"; else distType = "Delta"; end

% if distType == "Delta"
%     histColor = getColors().c_delta;
% else
%     histColor = getColors().c_uniform;
% end
distType = "";
histColor = [1 1 1];


index = struct ...
        ( minmaxrange_cs = find(unit.neuron(1).bins_cs >= -0.25 & unit.neuron(1).bins_cs <= 0.75) ...
        , maxmaxrange_cs = find(unit.neuron(1).bins_cs >= -0.25 & unit.neuron(1).bins_cs <= 0.6)...
        , base_idx_cs = find(unit.neuron(1).bins_cs >= -0.1 & unit.neuron(1).bins_cs <= -0.02) ...
        , complex_spike_range = find(unit.neuron(1).bins_cs >= 0 & unit.neuron(1).bins_cs <= 0.2)...  
        );
        
end