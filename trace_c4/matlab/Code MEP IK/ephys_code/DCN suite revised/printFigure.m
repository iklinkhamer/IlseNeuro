function printFigure(f, fileName, kwargs)
    arguments
        f (1,1) matlab.ui.Figure
        fileName (1,1) string
        kwargs.folder(1,1) string = pwd;
        kwargs.formats (1,:) string = ["png", "epsc"];
        kwargs.renderer (1,1) string = "painters";
        kwargs.size (1,2) {isnumeric}
    end
    
    function restoreFn = resize(fig, size)
        
        cachedUnits = fig.Units;
        cachedSize = fig.Position([3 4]);
        
        fig.Units = 'pixels';
        fig.Position([3 4]) = size;
        
        restoreFn = @() restoreSize(fig, cachedUnits, cachedSize);
        
        function restoreSize(fig, units, size)
            fig.Units = units;
            fig.Position([3 4]) = size;
        end
        
    end
    
    if isfield(kwargs, "size")
        cleanup = onCleanup(resize(f, kwargs.size));
    end
    
    arrayfun ...
        ( @(format) ...
            print ...
                ( f ...
                , sprintf("-d%s", format) ...
                , sprintf("-%s", kwargs.renderer) ...
                , fullfile ...
                    ( kwargs.folder ...
                    , fileName ...
                    ) ...
            ) ...
        , kwargs.formats ...
        );
    
end