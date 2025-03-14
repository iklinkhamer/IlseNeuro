mouse_names = ["Iowa", "Jackson", "Flint", "Greene", "Houston",  "Missouri", "Newark", "Orleans", "Pittsburg","Quimper", "Dallas", "Uppsala", "ReserveMouse3", "Seattle", "Madrid", "Zurich", "Yosemite", "Lisbon", "York", "Xanthi",  "Ana2", "Ana4", "Ana5", "Copenhagen", "Istanbul", "Kyiv", "Porto", "Queens", "Reno", "Rotterdam", "Tallinn", "Venice", "Willemstad", "Zachary", "Amsterdam"];
%mouse_names = ["Ana4"];
% done = []

for m = 1:length(mouse_names)
    try
        mouseName = mouse_names(m);
        inspectMouseUnits(mouseName, directory = fullfile(Env.getBayesLabUserRoot, "/TraceExperiments/AnalysisOutput/c4 results stats", mouseName), results_subfolder= "")
        JkUtils.memoize.clearCacheSystem();
    catch
        disp("Catch")
    end
end
