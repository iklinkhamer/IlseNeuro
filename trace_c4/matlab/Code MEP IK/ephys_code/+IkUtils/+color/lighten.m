%% I.K. 1-6-24
function color = lighten(color, multiplication_factor)
arguments
    color = [1 1 1];
    multiplication_factor = 1;
end
    color = color*0.7*multiplication_factor;            

end