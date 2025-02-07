%% IK 8-4-24
function stamp = nameDateStampFiles(kwargs)
arguments
    kwargs.mcode = '';
    kwargs.s = '';
    kwargs.file_pattern = '';
    kwargs.file_extension = '';
end
mcode = kwargs.mcode;
s = kwargs.s;
file_pattern = kwargs.file_pattern;
file_extension = kwargs.file_extension;
mice = defaultMice();

if length(split(mcode, '-')) < 3 && length(split(mcode, '.')) < 3 && ~ismember(mcode, ["Shank2KOMUT", "Shank2KOWT"])
    charmcode = char(mcode);
    midcode = string(charmcode(3:end));
    miceNames = [mice.code];
    miceSplitNames1 = split(miceNames(1:31), '-');
    miceMidNames1 = miceSplitNames1(:,:,2);
    miceMidNames1 = erase(miceMidNames1, 'MI');
    miceEndNames1 = miceSplitNames1(:,:,3);
    listMidCodes = miceMidNames1 + miceEndNames1;
    miceSplitNames2 = split(miceNames(32:end-2), '.');
    miceMidNames2 = miceSplitNames2(:,:,2);
    miceEndNames2 = miceSplitNames2(:,:,3);
    listMidCodes = [listMidCodes, miceMidNames2 + miceEndNames2];
    mcode = mice(listMidCodes == midcode).code;
end

mouse_idx = find([mice(:).code] == mcode);

datePattern = '^\d{6,7}$'; 
if (ischar(s) || isstring(s)) && ~isempty(regexp(s, datePattern, 'once')) && regexp(s, datePattern) %(strlength(s) == 6 || strlength(s) == 7)
    date = s;
elseif ~isempty(s) && ~isempty(mice(mouse_idx).ephysdates)
    date = mice(mouse_idx).ephysdates(s);
    if s == 3 && mice(mouse_idx).ephysdates(s) == mice(mouse_idx).ephysdates(s-1)
        date = date + "_2";
    end
elseif (isinteger(s) || isfloat(s)) && s <= IkUtils.getParams().n_sessions
    date = IkUtils.getParams().s(s);
else
    date = s;
end

if ~isempty(file_pattern)
    stamp = sprintf("%s", file_pattern);
end
if ~isempty(mcode)
    if exist("stamp", 'var')
        stamp = sprintf("%s_%s", stamp, mcode);
    else
        stamp = sprintf("%s", mcode);
    end
end
if ~isempty(date)
    if exist("stamp", 'var')
        stamp = sprintf("%s_%s", stamp, string(date));
    else
        stamp = sprintf("%s", string(date));
    end
end
if ~isempty(file_extension)
    if exist("stamp", 'var')
        stamp = sprintf("%s%s", stamp, file_extension);
    else
        stamp = sprintf("%s", file_extension);
    end
end
end