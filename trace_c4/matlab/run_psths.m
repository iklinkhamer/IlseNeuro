<<<<<<< HEAD
mouse_names = ["Quimper", "Uppsala", "ReserveMouse3", "Seattle", "Madrid", "Zurich", "Yosemite", "Lisbon", "York", "Xanthi", "Amsterdam", "Ana2", "Ana4", "Ana5", "Copenhagen", "Flint", "Greene", "Houston", "Iowa", "Istanbul", "Jackson", "Kyiv", "Missouri", "Newark", "Orleans", "Pittsburg", "Porto", "Queens", "Reno", "Rotterdam", "Tallinn", "Venice", "Willemstad", "Zachary", "Dallas"];
% mouse_names = ["Ana4"];
=======
mouse_names = [ "Jackson", "Flint", "Greene", "Houston", "Iowa", "Missouri", "Newark", "Orleans", "Pittsburg","Quimper", "Uppsala", "ReserveMouse3", "Seattle", "Madrid", "Zurich", "Yosemite", "Lisbon", "York", "Xanthi", "Amsterdam", "Ana2", "Ana4", "Ana5", "Copenhagen", "Istanbul", "Kyiv", "Porto", "Queens", "Reno", "Rotterdam", "Tallinn", "Venice", "Willemstad", "Zachary", "Dallas"];
%mouse_names = ["Ana4"];
>>>>>>> 17f49b4c80b9a6c2c8f54ff8dbb120a15e821c5b

for m = 1:length(mouse_names)
    try
        mouseName = mouse_names(m);
        inspectMouseUnits(mouseName)
    catch
        disp("Catch")
    end
end
