% TODO: Deprecate this function

function contamination = contaminationFilter(sessionsData)

    contamination = arrayfun ...
        ( @(session) computeContaminationForSession(session.neuron) ...
        , sessionsData ...
        );
%     contamination = struct;
%     
%     for s = 1: length(sessionsData)
%         unit = sessionsData(s);
%         
%         for idx = 1 : length(unit.neuron)
%             %% Finding the total contamination for a single neuron
%             
%             spike_intervals = diff(unit.neuron(idx).st) * 1000;
%             contamination(s).neuron(idx).spike_intervals = spike_intervals;
%             short_spike_intervals_mask = spike_intervals < 5;
%             
%             contamination(s).neuron(idx).percentage = sum...
%                 (short_spike_intervals_mask)/length(spike_intervals) * 100;
%         end
%     end
    
end

function contamination = computeContaminationForSession(neurons)
    
    contamination.neuron = arrayfun ...
        ( @computeContamination ...
        , neurons ...
        );
    
end

%% Julius: this commented section can probably be thrown out?

% histogramrange = find(unit.neuron(complexIdx).RasterXY_cs(1,:) >= -0.35 & unit.neuron(complexIdx).RasterXY_cs(1,:) <= 0.85);


%         raster_order_spikes = unit.neuron(complexIdx).RasterXY_cs(2,histogramrange);
%         raster_order_spikes_odd = raster_order_spikes(:,1:2:end);

%         spike_interval_counts = histcounts(spike_intervals,'BinWidth', 5, 'BinLimits', [0,1200],'FaceColor','b','EdgeColor','none','FaceColor',getColors().c_blue);
%         sum_spike_intervals = sum(spike_interval_counts.Values);
