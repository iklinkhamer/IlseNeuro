% 'fade' the current view of `ax` by drawing a partially transparent overlay
%
function [ax, overlay] = fadeAxes(ax, kwargs)
    arguments
        ax
        kwargs.color(1,3) {isnumeric} = [.6 .6 .6];
        kwargs.alpha(1,1) {isnumeric} = .3;
    end
    
    xRange = ax.XLim;
    x1 = xRange(1);
    x2 = xRange(2);
    yRange = ax.YLim;
    y1 = yRange(1);
    y2 = yRange(2);
    
    overlay = patch ...
        ( ax ...
        , [x1 x2 x2 x1], [y1 y1, y2 y2] ...
        , kwargs.color ...
        , EdgeColor = 'none' ...
        , FaceAlpha = kwargs.alpha ...
        );
    
end