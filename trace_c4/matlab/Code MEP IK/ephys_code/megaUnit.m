function megaUnit(kwargs)
arguments
    kwargs.loadDataFromDATA = true;
    kwargs.convertAllSessions = true;
    kwargs.convertAllMice = false;
    kwargs.startMouse = 1;
end

P = IkUtils.getParams();
% if kwargs.convertAllMice == 0
%
%     fprintf("\n-------------------------------------")
%     fprintf("\nWhich mouse would you like to curate?")
%     fprintf("\n-------------------------------------\n")
%
%     mouseCodes = arrayfun ...
%         ( @(mouse) string(mouse.code) ...
%         , defaultMice() ...
%         );
%     mouseNames = arrayfun ...
%         ( @(mouse) string(mouse.name) ...
%         , defaultMice() ...
%         );
%
%     [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
%     mcode = mouseCodes(mname_idx);
%     mousenames = mcode;
% else
%     mousenames = P.mouseList(kwargs.startMouse:end);
% end

sprintf("What group do you want to make a struct for?")
options = ["Shank2MUT", "Shank2WT"];
[group, group_idx] = IkUtils.do_prompt_select_option(options);

if group == "Shank2MUT"
%Shank2 MUT 
    mousenames = ["20-MI19442-03", "20-MI19442-05", "20-MI19442-08", "21-MI10159-01", "21-MI10159-02", "21-MI10532-07", "21-MI16091-04", "21-MI16091-03", "21-MI16183-05", "22-MI10447-09", "22-MI11756-06", "22-MI13134-01", "22-MI13989-07", "22-MI14020-07", "22-MI14020-08", "MI23.00102.03", "MI23.00244.04", "MI23.01412.04"];
elseif group == "Shank2WT"
    mousenames = ["20-MI19442-06", "21-MI10159-06", "21-MI10532-03", "21-MI16183-03", "22-MI10447-05", "22-MI10447-06", "22-MI11756-07", "22-MI13134-03", "22-MI13989-05", "22-MI14020-03", "MI23.00102.01", "MI23.01047.03", "MI23.01712.07"];
else
    disp("Group not recognized");
    return
end

if kwargs.convertAllSessions == 0
    prompt2 = "Enter session number: ";
    sessions = input(prompt2);
else
    sessions = 1:P.n_sessions;
end


n1 = 1;
n2 = 1;
n3 = 1;
for mname = mousenames
    structFilename = "StructEphysData.mat";
    try
        filename_try = fullfile(mname, P.s(1), structFilename);
        filename = which(filename_try);
        path_ = fileparts(filename);
        path = erase(path_, P.s(1));
    catch
        disp("Ephys Data Struct not found")
        data = [];
        continue
    end
    for s = sessions
        fprintf("Mouse: %s, %s  Session: %d \n", mname, P.mouseNames(find(P.mouseList == mname)), s);

       

        filename_try = fullfile(mname, P.s(s), structFilename);
        filename = which(filename_try);

        if exist(filename, 'file') == 2
            path_ = fileparts(filename);
            path = erase(path_, P.s(1));
        elseif s > 1
            continue
        else
            disp("Ephys Data Struct not found")
            %data = [];
            continue
        end


        load(filename);
        if length(unit.neuron) > 0
            unit.irc_trueID_all_neurons = unit.irc_trueID;
            unit = rmfield(unit, 'irc_trueID');
            for i = 1:length(unit.neuron)
                fields_struct = fields(unit);
                fields_struct = fields_struct(2:end);
                %                 for j = 1:length(fields_struct)
                %                     f = fields_struct{j};
                %                     unit.neuron(i).(f) = unit.(f);
                %                 end
                unit.neuron(i).mname = unit.mname;
                unit.neuron(i).name = unit.name;
                unit.neuron(i).session = unit.session;
                if s == 1
                    n = n1;
                elseif s == 2
                    n = n2;
                elseif s == 3
                    n = n3;
                end
                megaUnit(s).neuron(n) = unit.neuron(i);
                       if s == 1
                    n1 = n1 + 1;
                elseif s == 2
                    n2 = n2 + 1;
                elseif s == 3
                    n3 = n3 + 1;
                end

            end
        else
            fprintf("Data file not found for session %d\n", s)
            continue
        end

    end

end

disp("Saving...")
basepath = fileparts(fileparts(path));
savepath = fullfile(basepath, group);
for s = 1 : P.n_sessions
    try
        unit.neuron = megaUnit(s).neuron;

        stamp= nameDateStampFiles(mcode = group, file_pattern = 'StructEphysData', file_extension = '.mat');
        save(fullfile(savepath,P.s(s),stamp),'unit')
        fprintf("session %d saved", s)
    catch
        fprintf("session %d not saved\n", s)
    end
end

end