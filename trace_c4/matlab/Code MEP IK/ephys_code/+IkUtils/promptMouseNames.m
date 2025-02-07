function mouseCodes = promptMouseNames(mInput)

allMice = 0;
if nargin < 1
    fprintf("\n-------------------------------------")
    fprintf("\nWhich mouse would you like to curate?")
    fprintf("\n-------------------------------------\n")

    mouseCodes = arrayfun ...
        ( @(mouse) string(mouse.code) ...
        , defaultMice() ...
        );
    mouseNames = arrayfun ...
        ( @(mouse) string(mouse.name) ...
        , defaultMice() ...
        );

    mouseNames = mouseCodes + "     " + mouseNames;
    mouseNames(end+1) = "All mice";
    [mname, mname_idx] = IkUtils.do_prompt_select_option(mouseNames);
    try
        mouseCodes = mouseCodes(mname_idx);
    catch
        allMice = 1;
    end
    if allMice
        mouseCodes = arrayfun ...
            ( @(mouse) string(mouse.code) ...
            , defaultMice() ...
            );
    end
else
    mouseCodes = mInput;
end
end