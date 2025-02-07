function axs = initPlots(dim, fig)
arguments
    dim = [1,1];
    fig = figure("Position", [0 2000 1800 1000]);
end

fig;
for r = 1:dim(1)
    for c = 1:dim(2)

        axs(r,c) = subplot(dim(1), dim(2), (r-1)*dim(2)+c, Parent = fig);
        axs(r,c).TickDir = 'out';
        axs(r,c).NextPlot = 'add';
        axs(r,c).Box = false;
    end
end

end