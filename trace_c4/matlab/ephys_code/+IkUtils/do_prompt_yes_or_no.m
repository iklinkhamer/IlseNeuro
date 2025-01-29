function choice = do_prompt_yes_or_no(prompt)

fprintf("%s\n",prompt)
done = 0;
while ~done
prompt2 = "Type y for yes or n for no: ";

answer = string(input(prompt2, "s"));
if answer == 'y'
    choice = 1;
    done = 1;
elseif answer == 'n'
    choice = 0;
    done = 1;
else
    fprintf("\n No possible response detected. Type either y or n.");
end
end