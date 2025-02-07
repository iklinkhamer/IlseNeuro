function copyEphysDataFromNAS(mInput)
ephysSources = [...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-02_14-02-50_1/Record Node 101/"...
    "/run/user/1664601583/gvfs/afp-volume:host=Badura_NAS.local,user=i.klinkhamer,volume=Data2/Ongoing-Projects_2_(Backup-BaduraStore)/Eyeblink-EPHYS-ASD-Project/Ilse/Data/Eyeblink + Ephys/Ephys-data/13989-04M_2023-02-03_14-40-18_1/Record Node 101/"...
    ];
ephysDest = "/home/i.klinkhamer/Documents/Data/ephysData/";
for s = 1 : length(ephysSources)
    ephysSource = ephysSources(s);
    fileListEphys = dir(ephysSource);
    fileListEphys = fileListEphys(~[fileListEphys.isdir]);
    erasePart = fileparts(fileparts(fileparts(ephysSource)));
    ephysDestNew = fullfile(ephysDest, erase(ephysSource, erasePart));

    if ~exist(ephysDestNew, 'dir')
        mkdir(ephysDestNew);
    end

    for f = 1:length(fileListEphys)
        EphysSourceFile = fullfile(ephysSource, fileListEphys(f).name);
        EphysDestFile = fullfile(ephysDestNew, fileListEphys(f).name);
        if ~exist(EphysDestFile, 'file')
            copyfile(EphysSourceFile, EphysDestFile);
        end
    end
end
end