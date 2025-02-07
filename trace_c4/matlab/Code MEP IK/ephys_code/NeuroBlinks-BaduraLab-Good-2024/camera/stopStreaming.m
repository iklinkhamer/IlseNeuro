function stopStreaming(handles)

% Get rid of start/stop streaming and just run newFrameCallback all the time?

set(handles.togglebutton_stream,'Value',0);
set(handles.togglebutton_stream,'String','Start Streaming')

if isfield(handles, 'pwin')
    setappdata(handles.pwin(1),'UpdatePreviewWindowFcn',[]);
end