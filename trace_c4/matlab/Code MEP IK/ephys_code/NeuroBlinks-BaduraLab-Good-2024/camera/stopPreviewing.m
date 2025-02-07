function stopPreviewing(cameras)
% Pulled this out as a function so it can be called from elsewhere

gui = getappdata(0, 'gui'); 
handles = guidata(gui.maingui);

for i=1:length(cameras)
    stoppreview(cameras(i));
end

set(handles.pushbutton_StartStopPreview,'Value',0);
set(handles.pushbutton_StartStopPreview,'String','Start Preview')