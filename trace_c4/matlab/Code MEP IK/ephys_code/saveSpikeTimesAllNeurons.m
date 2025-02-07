%% I.K. 17/01/2024   Taking Spikesorting data outcome and putting it into a nice structure to use.
function saveSpikeTimesAllNeurons()

P = IkUtils.getParams();
app = getUImain(scriptname=mfilename);

while true
    uiwait(app.UIFigure)
    if ~isvalid(app)
        break; % If the UI figure has been closed, exit the loop
    end
    appData = evalin('base', 'appData');

    makeSpikeRasters = appData.makeSpikeRasters;
    mousecodes = appData.mousecodes;
    sessions = appData.sessions;

    for mcode = mousecodes
        no_neurons_count = 0;
        continue_mark = false;

        fprintf("Mouse: %s, %s\n", mcode, P.mouseNames(find(P.mouseList == mcode)))
        data = getData(mcode);
        if isempty(data)
            disp("Ephys data not found.")
            continue
        end
        modulationMasks = loadSimpleChannelModulationMasks(mcode);
        simpleMasks = loadSimpleMasks(mcode);
        trialspikesallsessions = [];
        if ~appData.singleSession
            sessions = 1 : length(modulationMasks);
        end

        for s = sessions
            fprintf("\nSession: %d \n", s);
            simpleMaskSession = simpleMasks{s};
            if sum(simpleMaskSession) == 0
                disp("No neurons in this session's simplemask.")
                continue
            end
            try
                if numel(data(s).neuron) == 0
                    disp("No neurons found in this session")
                    no_neurons_count = no_neurons_count + 1;
                    continue
                end
            catch
                continue_mark = true;
                continue
            end
            filepattern = "StructEphysData"; fileextension = ".mat";
            stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = filepattern, file_extension = fileextension);
            filename_try = fullfile(mcode, P.s(s), stamp);
            filename = which(filename_try);
            path = fileparts(filename);
            if isempty(path)
                disp("Ephys Data Struct not found")
                trialspikes = [];
                stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialspikes', file_extension = '.mat');
                save(fullfile(path, stamp),'trialspikes')
                continue
            end
            try
                if isempty(data(s).neuron)
                    disp("No neurons in this session.")
                    trialspikes = [];
                    stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialspikes', file_extension = '.mat');
                    save(fullfile(path,stamp),'trialspikes')
                    continue
                end
            catch
                continue
            end

            try
                sessionData = data(s).neuron(simpleMaskSession);
                n_neurons = length(sessionData);
            catch
                disp("Failed to determine length session data... Continuing to next session/mouse...")
                continue
            end

            trialspikes = struct;
            for n = 1 : n_neurons
                trial_spikes_cs_aligned = sessionData(n).trial_spikes_cs(:);
                trialspikes(n).spiketimes_cs_neuron = trial_spikes_cs_aligned;
                if length(trial_spikes_cs_aligned) < P.nCSplusRegTrials
                    for i = 1 : (P.nCSplusRegTrials-length(trial_spikes_cs_aligned))
                        trialspikes(n).spiketimes_cs_neuron = [];
                    end
                end
            end
            if ~isempty(trialspikes(1).spiketimes_cs_neuron)
                disp("saving session...")
                path = fullfile("/home/i.klinkhamer/Documents/Stuff/Folders", mcode, P.s(s));
                stamp = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'trialspikes', file_extension = '.mat');
                save(fullfile(path,stamp),'trialspikes')
                fprintf("saved\n")

                trialspikesallsessions = [trialspikesallsessions trialspikes];
            end            
            path2 = fileparts(path);
            if makeSpikeRasters
                n_neurons = length(trialspikes);
                n_trials = length(trialspikes(1).spiketimes_cs_neuron);
                bins = P.psthRanges.cs_full.min:P.spikeTimesBinW:P.psthRanges.cs_full.max;
                spikecounts = zeros(length(bins)-1, n_neurons*n_trials);
                i = 1;
                for n = 1 : n_neurons
                    spikes_neuron = trialspikes(n).spiketimes_cs_neuron;
                    if ~isempty(spikes_neuron)
                        for t = 1 : n_trials
                            spikes = trialspikes(n).spiketimes_cs_neuron(t).trial;
                            spikecounts(:,i) = histcounts(spikes,bins).';
                            i = i + 1;
                        end
                    end
                end
                if ~isempty(spikecounts)
                    disp("saving spike counts...")
                    stamp_scanat = nameDateStampFiles(mcode = mcode, s = s, file_pattern = 'spikecountsallneuronsalltrials', file_extension = '.mat');
                    save(fullfile(path,stamp_scanat),'spikecounts')
                    fprintf("saved\n")
                end
            end
        end
        if no_neurons_count == numel(sessions) || continue_mark
            continue
        end
        if ~isempty(trialspikesallsessions(1).spiketimes_cs_neuron)
            disp("saving all sessions...")
            stamp_scas = nameDateStampFiles(mcode = mcode, file_pattern = 'trialspikesallsessions', file_extension = '.mat');
            save(fullfile(path2,stamp_scas),'trialspikesallsessions')
            fprintf("saved\n")
        end

        if makeSpikeRasters
            n_neurons = length(trialspikesallsessions);
            n_trials = length(trialspikesallsessions(1).spiketimes_cs_neuron);
            bins = P.psthRanges.cs_full.min:P.spikeTimesBinW:P.psthRanges.cs_full.max;
            spikecounts = zeros(length(bins)-1, n_neurons*n_trials);
            i = 1;
            for n = 1 : n_neurons
                for t = 1 : n_trials
                    spikes = trialspikesallsessions(n).spiketimes_cs_neuron(t).trial;
                    spikecounts(:,i) = histcounts(spikes,bins).';
                    i = i + 1;
                end
            end
            if ~isempty(spikecounts)
                disp("saving all neurons, all trails, all sessions...")
                stamp_scanalas = nameDateStampFiles(mcode = mcode, file_pattern = 'spikecountsallneuronsalltrialsallsessions', file_extension = '.mat');
                save(fullfile(path2,stamp_scanalas),'spikecounts')
                fprintf("saved\n")
            end
        end

    end

    fprintf("\nLoop done, waiting for new UI input confirmation.\n")
end
end
