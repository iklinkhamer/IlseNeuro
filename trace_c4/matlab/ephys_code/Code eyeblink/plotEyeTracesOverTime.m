%% IK 8-4-24
function plotEyeTracesOverTime(mouseData, ID, fig, subcoords, day)
P = IkUtils.getParams();
lastDay = min(IkUtils.getParams().lastDay, sum(any(~isnan(mouseData.onsets), 2)));

if nargin < nargin('plotEyeTracesOverTime')
    day = lastDay;
end

axes(fig);
% set(gcf, 'Renderer', 'zbuffer'); % Change renderer to zbuffer
set(gcf, 'Renderer', 'opengl'); % Change renderer to opengl
subplot(subcoords(1), subcoords(2), subcoords(3))
hold on

mouseDatameanWF=mouseData.mean;
mouseDatameanWF(:,1)=NaN; %put first and last timepoint on NaN to prevent "curtain effect" at edges
mouseDatameanWF(:,P.nTimeSteps)=NaN;

colormap(winter)
for i=1:P.lastDay
    color(i,1:P.nTimeSteps)=i*(length(mouseDatameanWF)); %give each line a different color (instead of gradient for heigth profile)
end
try
    h1=waterfall(mouseData.ts,1:P.lastDay,mouseDatameanWF(1:day,:),color);
catch
    h1=waterfall(mouseData.ts,1:day,mouseDatameanWF(1:day,:),color(1:day, :));  
end

h1.LineWidth=2;
h1.FaceAlpha=0; %transparent faces


view([-0.5 -1.5 1]) %determines view angle
yticks(1:day);
xlim([-200 800]);

if length(ID) > 1
    Gene = ID(1);
    loc = ID(2);
    Type = ID(3);
    title(sprintf('Average traces day 1-%d: %s %s %s', day, Gene, loc, Type))
else
    title(sprintf('Average traces day 1-%d: %s', day, ID))
end

xlabel('Time from CS (ms)')
zlabel('fraction eyelid closure')
ylabel_handle = ylabel('training day');
set(ylabel_handle, 'Rotation', -45, 'HorizontalAlignment', 'right');
set(gca, "TickDir", "out")
box("off")
end