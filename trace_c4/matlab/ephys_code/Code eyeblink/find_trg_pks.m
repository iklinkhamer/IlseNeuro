% Find the rising phase peaks of recorded trigger data.
% (c) Si-yang Yu @ PSI 2018

% Function input:  trg,   the raw recording of trigger data
%                  v_trg, trigger voltage used (V)
%                  t,     recording timestamps, must be the same length as trg
%                  dur,   experiment trial duration (ms)
% Function output: locs,  trigger rising phase peak times


function locs = find_trg_pks(trg,v_trg,t,dur)
    trg = round(trg);
    v_trg = round(v_trg) / 2;
    dur = dur / 1000;
    [~,locs] = findpeaks(trg,t,'MinPeakHeight',v_trg,'MinPeakDistance',dur);
end