%% I.K. 1-6-24
clear;
mouseCodes = arrayfun ...         
    ( @(mouse) string(mouse.code) ...
    , defaultMice() ...
    );
mouseNames = arrayfun ...
    ( @(mouse) string(mouse.name) ...
    , defaultMice() ...
    );

ephysMiceMask = ~cellfun(@isempty,{defaultMice().ephysdates});
ephysMiceMask(end-1:end) = true;
mouseNames = mouseNames(ephysMiceMask);
mouseCodes = mouseCodes(ephysMiceMask);
mouseList = mouseCodes;
mouseList(end + 1) = "Shank2KOMUT";
mouseList(end + 1) = "Shank2KOWT";

days_ = 1:31;
for i = 1:length(days_)
    day = days_(i);
    if day < 10
        day = sprintf('0%d',day);
    else
        day = int2str(day);
    end
    days{i} = day;
end
months_ = 1:12;
for i = 1:length(months_)
    month = months_(i);
    if month < 10
        month = sprintf('0%d',month);
    else
        month = int2str(month);
    end
    months{i} = month;
end
years = 20:24;

defaultMiceList = defaultMice();


for i = 1:length(mouseList)
    numDates = length(defaultMiceList([defaultMiceList.code] == mouseList(i)).ephysdates);

    status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s", mouseList(i)));

    if numDates > 0
        status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S1", mouseList(i)));
    end
    if numDates > 1 || mouseList(i) == "Shank2MUT" || mouseList(i) == "Shank2KO"
        status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S2", mouseList(i)));
    end
    if numDates > 2 || mouseList(i) == "Shank2MUT" || mouseList(i) == "Shank2KO"
        status = mkdir(sprintf("/home/i.klinkhamer/Documents/Stuff/Folders/%s/S3", mouseList(i)));
    end
end