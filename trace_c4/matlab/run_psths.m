mouse_names = [ "Jackson", "Flint", "Greene", "Houston", "Iowa", "Missouri", "Newark", "Orleans", "Pittsburg","Quimper", "Uppsala", "ReserveMouse3", "Seattle", "Madrid", "Zurich", "Yosemite", "Lisbon", "York", "Xanthi", "Amsterdam", "Ana2", "Ana4", "Ana5", "Copenhagen", "Istanbul", "Kyiv", "Porto", "Queens", "Reno", "Rotterdam", "Tallinn", "Venice", "Willemstad", "Zachary", "Dallas"];
%mouse_names = ["Ana4"];

for m = 1:length(mouse_names)
    try
        mouseName = mouse_names(m);
        inspectMouseUnits(mouseName)
    catch
        disp("Catch")
    end
end