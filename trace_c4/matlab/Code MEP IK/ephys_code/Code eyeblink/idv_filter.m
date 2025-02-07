% Apply different combination of LP/HP/NOTCH filters for Intan recordings.
% (c) Si-yang Yu @ PSI 2018

% Function input:  v_in,        defined as an array of voltage data
%                  apply_lp,    applying lowpass filter if this value NOT EQUAL to 0
%                  lp_freq,     lowpass filter frequency
%                  apply_hp,    applying highpass filter if this value NOT EQUAL to 0
%                  hp_freq,     highpass filter frequency
%                  apply_notch, applying notch filter if this value NOT EQUAL to 0
%                  notch_freq,  notch filter frequency
% Function output: v_out,       output array of filtered voltage data

function v_out = idv_filter(v_in,apply_lp,lp_freq,apply_hp,hp_freq,apply_notch,notch_freq)
    fs = evalin('caller','frequency_parameters.amplifier_sample_rate');
  % Converting inputs
    apply_lp = logical(apply_lp);
    apply_hp = logical(apply_hp);
    apply_notch = logical(apply_notch);
  % Verify inputs
    if apply_lp && apply_hp
        if lp_freq <= hp_freq
            error('Lowpass Filter freqency MUST be greater than Highpass Filter freqency!');
        end
    end
    if apply_lp && apply_notch
        if lp_freq <= notch_freq
            warning('Notch Filter freqency greater than Lowpass Filter freqency, Notch MAY NOT effective!');
        end
    end
    if apply_hp && apply_notch
        if notch_freq <= hp_freq
            warning('Notch Filter freqency lesser than Highpass Filter freqency, Notch MAY NOT effective!');
        end
    end
  % Processing all filters
    if apply_lp
        fNorm_low = lp_freq / (fs / 2);
        [b_low,a_low] = butter(10, fNorm_low, 'low');
        v_out = filtfilt(b_low, a_low, v_in);
    end
    if apply_hp
        fNorm_high = hp_freq / (fs / 2);
        [b_high,a_high] = butter(3, fNorm_high, 'high');
        v_out = filtfilt(b_high, a_high, v_in);
    end
    if apply_notch   
        fNorm_notch = notch_freq / (fs / 2);
        notch_bandw = fNorm_notch / 10;
        [b_notch, a_notch] = iirnotch(fNorm_notch, notch_bandw);
        v_out = filtfilt(b_notch, a_notch, v_in);
    end
end