% Callback for camera2
function newFrameCallbackCam2(obj, event, himage)

    persistent last_pupil_diameter
    persistent last_pupil_center
    
    gui = getappdata(0, 'gui'); 
    handles = guidata(gui.maingui);
    % cameras = getappdata(0, 'cameras');
    metadata = getappdata(0, 'metadata');  
    
    
    % --- eye trace ---
    wholeframe = event.Data;
    roi = wholeframe .* uint8(metadata.cam(2).mask);
    eyelidpos = sum(roi(:) >= 255 * metadata.cam(2).thresh);
    
    % Pupil
    pupil_thresh = str2double(handles.edit_pupil_thresh.String);
    % Value of zero means don't stream it; otherwise, extract pupil diameter
    if pupil_thresh > 0
        % Extract pupil diameter
        % Use eyelid ellipse ROI to guess pupil center
        pos = metadata.cam(2).winpos;
        ellipse_center = [pos(1) + pos(3)./2, pos(2) + pos(4)./2];
        % Need to add ability to fine tune parameters to this function
        roi_x = round(pos(3) ./ 2);
        roi_y = round(pos(4) ./ 2);
    %     max_radius = round(min(roi_x.*2, roi_y.*2));
    %     max_radius = round(min(roi_x, roi_y));
        max_radius = 35;
        [area, center] = getPupilArea(wholeframe, ellipse_center, [roi_x, roi_y], [max(1, round(max_radius ./ 3)), max_radius], pupil_thresh);
        r = sqrt(area ./ pi);
        pupil_diameter = 2 .* r;
    
        % Add circle to frame ('CData')
        if ~isnan(r)
            % Just reuse last value if it changed too fast (due to blink or whiskers getting in ROI)
            % if abs(pupil_diameter-last_pupil_diameter) > 5
            %     pupil_diameter = last_pupil_diameter;
            %     center = last_pupil_center;
            %     r = last_pupil_diameter./2;
            % end
    
            circle_color = 'blue';
            last_pupil_diameter = pupil_diameter;
            last_pupil_center = center;
        else
            center = last_pupil_center;
            r = last_pupil_diameter ./ 2;
            pupil_diameter = last_pupil_diameter;
            circle_color = 'red';
        end
        
        circle_params = [center, r];
            
        F = insertShape(event.Data, 'Circle', circle_params, 'Color', circle_color, 'LineWidth', 2);
    
        % Save pupil diameter to disk
        % NOT IMPLEMENTED YET

    else
        F = event.Data;
    end
