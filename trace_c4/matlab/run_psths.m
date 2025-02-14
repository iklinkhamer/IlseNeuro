mouse_names = ["Quimper", "Uppsala", "ReserveMouse3", "Seattle", "Madrid", "Zurich", "Dallas", "Yosemite", "Lisbon", "York", "Xanthi", "Amsterdam", "Ana2", "Ana4", "Ana5", "Copenhagen", "Flint", "Greene", "Houston", "Iowa", "Istanbul", "Jackson", "Kyiv", "Missouri", "Newark", "Orleans", "Pittsburg", "Porto", "Queens", "Reno", "Rotterdam", "Tallinn", "Venice", "Willemstad", "Zachary"];


for m = 1:length(mouse_names)
    try
        mouseName = mouse_names(m);
        inspectMouseUnits(mouseName)
    catch
        disp("Catch")
    end
end