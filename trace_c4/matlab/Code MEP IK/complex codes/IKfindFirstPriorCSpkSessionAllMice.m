function [indices, names] = IKfindFirstPriorCSpkSessionAllMice()
    
    mice = IkUtils.filter ...
        ( @(m) any(strcmpi(m.type, ["uniform", "deun"])) ...
        , ilseFigsMice() ...
        );
    
    names = [mice.name];
    indices = arrayfun ...
        ( @findFirstPriorCSpkSession ...
        , names ...
        );
    
    figure;
    box off
    set(0, 'defaultFigureRenderer', 'painters')
    set(gca,'TickDir','out')
    counts = ...
        [ 50.6  87.4
        27.5    44.7
        42  47.2
            ];
        
    bar(counts);
%     counts = ...
%         [ 2 
%          3 
%          1 
%         ];
%  x = [1,1,1];
%  y = [1,2,0];
% %     
% %     ax = IkUtils.initPlots();
% %     [y, b] = histogram(ax, [indices(1:3).', indices(4:6).']);
% %     bar(b,y, 'grouped')
%     %histogram(ax, indices(4:6))
%     figure(2);
% %     histogram([x.';y.'])
%     bar(counts)
%     title("Indices of first uniform session with prior complex spikes")
     
    
end