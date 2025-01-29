function togglePreview(handles)

cameras = getappdata(0, 'cameras');
metadata = getappdata(0, 'metadata');
gui = getappdata(0, 'gui');

% if ~isfield(metadata.cam(1),'fullsize')
%     metadata.cam(1).fullsize = [0 0 640 480];
% end
% metadata.cam(1).ROIposition=camera.ROIposition;

% Start/Stop Camera
if strcmp(get(handles.pushbutton_StartStopPreview, 'String'), 'Start Preview')
    % Camera is off. Change button string and start camera.
    set(handles.pushbutton_StartStopPreview, 'String', 'Stop Preview')
    % Send camera preview to GUI
    for i = 1:length(cameras)

        % Need to change this code so we don't have to delete and re add roi ellipse and trialtimecounter in order to put them on top of newly created pwin image

        % This line here is probably responsible for causing the "invalid or deleted object" error during newFrameCallback
        % delete(gui.cameraAx(i).Children)
        if isfield(handles,'pwin')
            if length(handles.pwin) == length(cameras)  % Only delete it if it has already been created
                delete(handles.pwin(i))
            end
        end

        imx = metadata.cam(i).ROIposition(1) + [1:metadata.cam(i).ROIposition(3)];
        imy = metadata.cam(i).ROIposition(2) + [1:metadata.cam(i).ROIposition(4)];
%         imx = 0 + [1:metadata.cam(i).ROIposition(3)];
%         imy = 0 + [1:metadata.cam(i).ROIposition(4)];
        handles.pwin(i) = image(imx, imy, zeros(metadata.cam(i).ROIposition([4 3])), 'Parent', gui.cameraAx(i));
        
        % Turn off warnings to suppress message about trigger mode
        warning('off')
        % -----------------------------------------------------------------
        % This is a temporary workaround for what appears to be a bug with
        % the Ace camera, where the ROIPosition gets reset upon starting
        % the preview
        oldROIPosition = cameras(i).ROIPosition;
        % -----------------------------------------------------------------
        
        preview(cameras(i), handles.pwin(i));
        
        % -----------------------------------------------------------------
        % This is a temporary workaround for what appears to be a bug with
        % the Ace camera, where the ROIPosition gets reset upon starting
        % the preview
        cameras(i).ROIPosition = oldROIPosition;
        % -----------------------------------------------------------------
        
        warning('on')
        % set(gui.cameraAx(i), 'XLim', 0.5 + metadata.cam(i).ROIposition([1 3]))
        % set(gui.cameraAx(i), 'YLim', 0.5 + metadata.cam(i).ROIposition([2 4]))
        set(gui.cameraAx(i), 'XLim', [imx(1) imx(end)])
        set(gui.cameraAx(i), 'YLim', [imy(1) imy(end)])
        hp = findobj(gui.cameraAx(i),'Tag','roipatch');  
        delete(hp)
%         if ~isempty(hp)
%             uistack(hp, 'top')
%         end

        % if isfield(handles,'XY')
        %     handles.roipatch = patch(handles.XY(:, 1), handles.XY(:, 2), 'g', 'FaceColor', 'none', 'EdgeColor', 'g', 'Tag', 'roipatch');
        % end
        
    end

    ht = findobj(gui.cameraAx(1), 'Tag', 'trialtimecounter');
    delete(ht)
    
    % Put countdown on first preview window
    % axes(gui.cameraAx(1))
    ht = findobj(gui.cameraAx(1),'Tag','trialtimecounter'); 
    if isempty(ht)
        handles.trialtimecounter = text(gui.cameraAx(1), handles.cameraAx1.XLim(2)-10, handles.cameraAx1.YLim(2)-10, ...
            ' ', 'Color', 'c', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom', ...
            'Visible', 'Off', 'Tag', 'trialtimecounter', 'FontSize', 18);
    end

else
    % Camera is on. Stop camera and change button string.
    set(handles.pushbutton_StartStopPreview,'String','Start Preview')
    stopPreviewing(cameras);
end

setappdata(0,'metadata', metadata);
guidata(gui.maingui, handles)

drawnow         % Seems necessary to update appdata before returning to calling function