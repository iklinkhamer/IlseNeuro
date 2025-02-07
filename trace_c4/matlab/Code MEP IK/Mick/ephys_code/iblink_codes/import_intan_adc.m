% Reads an on board ADC data file and creates a voltage (V) struct array.
% (c) Si-yang Yu @ PSI 2018

% Function input:  adc_ch,  defined as the used on board ADC channel(s) for recording, read all existing if empty
% Function output: adc_sgn, defined as an array of struct (2 fields, defined below)
%                             nr: the number of used ADC channel
%                             v:  the votage signal of ADC channel (V)

function adc_sgn = import_intan_adc(adc_ch)
    warning('off','backtrace');
  % Proofreading of input
    adc_ch = round(adc_ch);
    adc_ch(adc_ch > 7) = 7;
    adc_ch(adc_ch < 0) = 0;
    adc_ch = unique(adc_ch);
  % Get related workspace data from caller
    try
        chip_channel = evalin('caller','[board_adc_channels.chip_channel]');
    catch
        chip_channel = [];
    end
  % Checking if input channel(s) in the recording
    if ~isempty(adc_ch)
        n = 1;
        while n <= size(adc_ch,2)
            if isempty(find(chip_channel == adc_ch(n),1))
                warning(strcat('On board ADC channel No.',num2str(adc_ch(n),'%02d'),' not recorded!'));
            end
            n = n + 1;
        end
        adc_ch = intersect(adc_ch,chip_channel);
    else
        disp('Empty ADC channel input! Reading all exsisting ADC channel(s)!');
        adc_ch = chip_channel;
    end
  % Execution session
    if ~isempty(adc_ch)
        n_adc_ch = size(adc_ch,2);
        i = 1;
        adc_sgn(n_adc_ch) = struct('nr',[],'v',[]);
        while i <= n_adc_ch
            file_no = num2str(adc_ch(i),'%02d');
            file_name = strcat('board-ADC-',file_no,'.dat');
            if (exist(file_name,'file'))
              % Reads a board ADC data file and creates an electrode voltage vector (V)
                fileinfo = dir(file_name);
                num_samples = fileinfo.bytes/2;
                fid = fopen(file_name,'r');
                v_adc = fread(fid,num_samples,'uint16') * 0.000050354;  % value_type = 'uint16'; mutiplier = 0.000050354 (V)
                fclose(fid);
                adc_sgn(i) = struct('nr',adc_ch(i),'v',v_adc);
            else
              % Assign an empty value to missing / unused ADC channel(s)
                adc_sgn(i) = struct('nr',adc_ch(i),'v',[]);
                warning(strcat('On board ADC channel No.',num2str(adc_ch(i),'%02d'),' file missing!'));
            end
            i = i + 1;
        end
    else
      % No ADC channel(s) defined or found
        adc_sgn = struct('nr',{},'v',{});
        warning('No existing on board ADC channel found in header file!');
    end
    warning('on','backtrace');
end