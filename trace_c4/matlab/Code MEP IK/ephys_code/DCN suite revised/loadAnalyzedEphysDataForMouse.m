%% I.K. 1-6-24
function sessions = loadAnalyzedEphysDataForMouse(mname, kwargs)
    arguments
        mname(1,1) string
        kwargs.writeToCache (1,1) logical = false;
        kwargs.readFromCache(1,1) logical = false;
    end
    
    %% MEMOIZE RETRIEVE
    if kwargs.readFromCache
        [sessions, ok] = IkUtils.memoize.fetch(mname);
        if ok
            return
        end
    end
    %%
%     
%     fileLists = getEphysFilesForMouse(mname);
%     
%     analyzedEphysFiles = fileLists.analyzedEphys;
    
    %sessions = getAnalyzedEphysData(analyzedEphysFiles);  % IK change
    sessions = getData(mname);
    
    %% MEMOIZE STORE
    if kwargs.writeToCache
        IkUtils.memoize.store(mname, sessions);
    end
    %%
    
end