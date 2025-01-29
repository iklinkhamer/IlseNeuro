function type = sessionType(sessionData)

%     nConditions = 2; %size(sessionData.neuron(1).AveFiring_cs,2);
defaultMiceList = defaultMice();    
type = defaultMiceList([defaultMice().code] == sessionData.mname).type;

%     switch nConditions
%         case 2
%             type = "Delta";
%         case 6
%             type = "Uniform";
%         otherwise
%             error ...
%                 ( "Could not assign distribution type: unknown number of conditions: %d" ...
%                 , nConditions ...
%                 )
%     end    

end