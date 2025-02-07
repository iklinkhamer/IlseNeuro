function [choice, choice_idx] = do_prompt_select_option(choices)
for i = 1:length(choices)
    fprintf("%d. %s\n", i, choices(i))
end
choice_made = false;
while choice_made == false

    choice_idx = input("", 's');
    choice_idx = str2num(choice_idx);
    if isempty(choice_idx) || choice_idx > length(choices)
        disp("Input not found in choice list. Please enter a number corresponding to one of the choices.")
    else
        choice = choices(choice_idx);
        break
    end
end
end
