% Given a mouse name, load its sessions from disk and return the indices
% corresponding to Shank2WT and Shank2MUT sessions as follows:
%   - Shank2WT sessions are those sessions that have type "Shank2WT" and that
%   occur before any Shank2MUT sessions have occured.
%   - Shank2MUT sessions are those sessions that form a contiguous list
%   of type "Shank2MUT" occuring after the first occurrence of type
%   "Shank2MUT".
%   - any sessions that do not follow either of the above criteria are
%   discarded.
%   - Both the Shank2WT indices and the Shank2MUT indices are potentially empty

% Copyright 2023 NarainLab
% Erasmus Medical Center, Rotterdam, The Netherlands
% v2023-01-21
% 
function sessionIdcs = splitSessionsByType(mname)    

    sessions = loadAnalyzedEphysDataForMouse(mname);
    nSessions = numel(sessions);
    sessTypes = arrayfun(@sessionType, sessions);
    
    firstShank2MUT = find(strcmpi("Shank2MUT", sessTypes), 1);
    if isempty(firstShank2MUT)
        Shank2MUTIdcs = [];
        Shank2WTIdcs = 1:nSessions;
    else
        Shank2WTAfterShank2MUTOffset = find(strcmpi("Shank2WT", sessTypes(firstShank2MUT+1:end)));
        if isempty(Shank2WTAfterShank2MUTOffset)
            lastShank2MUT = nSessions;
        else
            lastShank2MUT = firstShank2MUT + Shank2WTAfterShank2MUTOffset - 1;
        end
        Shank2MUTIdcs = firstShank2MUT : lastShank2MUT;
        Shank2WTIdcs = 1:(firstShank2MUT-1);
    end
    
%     type = getMouseType(mname);
    
%     switch lower(type)
%         case "Shank2WT"
%             sessionIdcs = struct ...
%                 ( Shank2WT = Shank2WTIdcs ...
%                 , Shank2MUT = [] ...
%                 );
%             
%         case "Shank2MUT"
%             sessionIdcs = struct ...
%                 ( Shank2WT = [] ...
%                 , Shank2MUT = Shank2MUTIdcs ...
%                 );
% 
%         case "deun"
            sessionIdcs = struct ...
                ( Shank2WT = Shank2WTIdcs ...
                , Shank2MUT = Shank2MUTIdcs ...
                );
%     end

    
%     if numel([Shank2WTIdcs Shank2MUTIdcs]) < numel(sessions)
%         keyboard
%     end
    
end