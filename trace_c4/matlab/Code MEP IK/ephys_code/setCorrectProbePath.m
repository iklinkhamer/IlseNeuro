function setCorrectProbePath(mname, s)
P = getParams();
probePath = P.prbPath;
path = fullfile(P.dirHome, mname, P.s(s));

if isfile(fullfile(path,"Config_h32_oe_bin.prm"))
    filename = fullfile(path, "Config_h32_oe_bin.prm");
elseif isfile(fullfile(path, "Config_h32_oe_dat.prm"))
    filename = fullfile(path, "Config_h32_oe_dat.prm");
elseif isfile(fullfile(path, "Config_h32_oe.prm"))
    filename = fullfile(path, "Config_h32_oe.prm");
else
    disp("Params file not found")
    return
end


setParamDataPath(filename, probePath);
end

function setParamDataPath(filename, probePath)
file  = fopen(filename,'r');
try
    file_contents = fread(file, 'char*1');
catch
    file_contents = fileread(filename);
end
file_contents_new = strrep(file_contents, "'DBC_3.1-64-H2_IK.prb'", sprintf("'%s'", probePath));
file = fopen(filename,'w');
fprintf(file, '%s', file_contents_new);
fclose(file);
end