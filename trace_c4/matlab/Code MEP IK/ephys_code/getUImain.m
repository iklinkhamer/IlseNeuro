%% IK 24-2-24
function app = getUImain(kwargs)
arguments
    kwargs.field = ''
    kwargs.value = ''
    kwargs.items = ''
    kwargs.scriptname = ''
end
app = uimain();
try
    mouseNames = [[defaultMice().name], "N/A"];
    mouseCodes = [[defaultMice().code], "N/A"];
catch
    IkUtils.addPaths();
    mouseNames = [[defaultMice().name], "N/A"];
    mouseCodes = [[defaultMice().code], "N/A"];
end
app.MouseDropDown.Items = mouseNames;
app.CodeDropDown.Items = mouseCodes;
app.StartmousenameDropDown.Items = mouseNames;
app.StartmousecodeDropDown.Items = mouseCodes;
if ~isempty(kwargs.field)
    for f = 1 : length(kwargs.field)
        field = kwargs.field(f);
        if ~isempty(kwargs.value)
            value = kwargs.value(f);
            app.(field).Value = value;
        end
        if ~isempty(kwargs.items)
            items = kwargs.items(f);
            app.(field).Items = items;
        end
    end
end
if ~isempty(kwargs.scriptname)
    appFields = fields(app);
    idx = 1;
    for c = 1 : numel(app.FunctionsettingsPanel.Children)
        if contains(class(app.FunctionsettingsPanel.Children(c)), "Panel")
            titles(idx) = string(app.FunctionsettingsPanel.Children(c).Title);
            idx = idx + 1;
        end
    end
    children = [];
    for f = 1 : numel(appFields)
        field = string(appFields(f));
        try
            if app.(field).Parent.Title == "Function settings"
                children = [children field];
            end
        catch
        end
    end
    for c = 1 : numel(children)
        child = children(c);
        if child == kwargs.scriptname
            app.(child).Visible = true;
        else
            app.(child).Visible = false;
        end
    end
else
    appFields = fields(app);
    for f = 1 : numel(appFields)
        field = appFields(f);
        try
            app.(field).Visible = true;
        catch
        end
    end
end
end