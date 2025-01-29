classdef UiState
    
    properties
        channelIdx
        
        masks
        
        reachedBottom
        
        reachedTop
        
        finishSession
    end
    
    methods
        function state = UiState(kwargs)
            arguments
                kwargs.channelIdx(1,1) {isnumeric}
                kwargs.masks(1,:) struct
                kwargs.exemplaryMask(1,:) {islogical}
                kwargs.reachedBottom(1,1) {islogical}
                kwargs.reachedTop(1,1) {islogical}
                kwargs.finishSession (1,1) logical
            end
            
            fields = string(fieldnames(kwargs));
%             values = struct2cell(kwargs);
            
            for field = fields(:)'
                state.(field) = kwargs.(field);
            end
%             cellfun ...
%                 ( @(field, value) set(state, field, value) ...
%                 , fields ...
%                 , values ...
%                 )
              
            
        end
    end
    
end