function ok = stopSession(handles)

cameras=getappdata(0,'cameras');
microController=getappdata(0,'microController');
timestamp = getappdata(0, 'timestamp');
metadata = getappdata(0, 'metadata');

pause(1)
% while (cameras(2).FramesAcquired ~= cameras(2).DiskLoggerFrameCount) 
%      pause(.1)
% end

timestamp.trial_time(timestamp.counter) = etime(clock, datevec(metadata.ts(1)));
trialTimes_cam2 = timestamp.trial_time - timestamp.trial_time(1);
ok = 0;

stopStreaming(handles)
stopPreviewing(cameras);

%FramesAcquired = cameras(2).FramesAcquired
%DiskLoggerFrameCount = cameras(2).DiskLoggerFrameCount 
name = sprintf('%s\\%s_%s.mat',metadata.folder, metadata.basename, 'timestampsCam2')
    save(name, 'trialTimes_cam2');
	
name2 = sprintf('%s\\%s_%s.mat',metadata.folder, metadata.basename, 'timestampsTotalCam2')
    save(name2, 'timestamp');
	
try
    fclose(microController);
    delete(microController);
    delete(cameras);
    rmappdata(0,'cameras');

    ok = 1;

catch err
    warning(err.identifier,'Problem cleaning up objects. You may need to do it manually.')
end
