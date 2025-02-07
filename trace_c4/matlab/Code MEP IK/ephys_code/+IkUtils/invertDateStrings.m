% Example date string in 'DD-MM-YYYY' format
mouseList = defaultMice();
dates = [mouseList.ephysdates];
% dateString = '15-03-2024';

for dateString = dates
    % Split the date string into day, month, and year
    parts = split(dateString, '-');

    try
        % Rearrange the parts to 'YYYY-MM-DD' format
        newDateString = sprintf('%s-%s-%s', parts{3}, parts{2}, parts{1});
    catch
        continue
    end

    % Display the converted date string
    disp(['Original date string: ', dateString]);
    disp(['Converted date string: ', newDateString]);
end