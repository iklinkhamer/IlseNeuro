%% I.K. 1-6-24
function p = Params()
p.fontName = 'Ariel';
p.lineWidth = 10;
p.figureForegroundColor = [0.5 0.5 0.5];

p.black = [
    0.00, 0.00, 0.00;  % Black
    0.10, 0.10, 0.10;  % Dark grey
    0.20, 0.20, 0.20;  % Grey
    0.30, 0.30, 0.30;  % Light grey
    0.40, 0.40, 0.40;  % Light grey
    0.50, 0.50, 0.50;  % Light grey
    0.60, 0.60, 0.60;  % Light grey
    0.70, 0.70, 0.70;  % Light grey
    0.80, 0.80, 0.80;  % Light grey
];

p.green = [
    0, 204, 0;
    0, 230, 0;
    26, 255, 26;
    51, 255,51;
    77, 255,77;
    102, 255,102;
    128, 255,128;
    153, 255,153;
    179, 255,179
    ] / 255;

% Define the color gradients
p.orange = [
    255, 107, 0;    % Dark Orange
    255, 123, 26;
    255, 139, 51;
    255, 155, 77;
    255, 171, 102;
    255, 186, 128;
    255, 201, 153;
    255, 217, 179;
    255, 232, 204    % Light Orange
] / 255;

p.red = [...
    204, 0, 0;
    230, 0, 0;
    255, 26, 26;
    255, 51, 51;
    255, 77, 77;
    255, 102, 102;
    255, 128, 128;
    255, 153, 153;
    255, 179, 179;
    ] / 255;



% Number of colors to interpolate
num_colors = size(p.orange, 1);

% Interpolate between the two color schemes
p.orange_red = zeros(num_colors, 3);

for i = 1:num_colors
    p.orange_red(i, :) = (p.orange(i, :) + p.red(i, :)) / 2;
end


% Define the color gradients
p.purple = [
    102, 0, 153;    % Dark Purple
    122, 28, 171;
    142, 56, 189;
    162, 84, 207;
    182, 112, 225;
    202, 140, 243;
    217, 163, 255;
    227, 183, 255;
    237, 203, 255    % Light Purple
]/255;

p.blue = [
    0, 0, 204;
    0, 0, 230;
    26, 26, 255;
    51, 51, 255;
    77, 77, 255;
    102, 102, 255;
    128, 128, 255;
    153, 153, 255;
    179, 179, 255
] / 255;



% Number of colors to interpolate
num_colors = size(p.purple, 1);

% Interpolate between the two color schemes
p.purple_blue = zeros(num_colors, 3);

for i = 1:num_colors
    p.purple_blue(i, :) = (p.purple(i, :) + p.blue(i, :)) / 2;
end


% Interpolate between the two color schemes
p.dark_red = zeros(num_colors, 3);

for i = 1:num_colors
    p.dark_red(i, :) = (p.red(i, :) + p.black(i, :)) / 2;
end

% Interpolate between the two color schemes
p.darker_red = zeros(num_colors, 3);

for i = 1:num_colors
    p.darker_red(i, :) = (p.dark_red(i, :) + p.black(i, :)) / 2;
end

% Interpolate between the two color schemes
p.dark_green = zeros(num_colors, 3);

for i = 1:num_colors
    p.dark_green(i, :) = (p.green(i, :) + p.black(i, :)) / 2;
end

% Interpolate between the two color schemes
p.darker_green = zeros(num_colors, 3);

for i = 1:num_colors
    p.darker_green(i, :) = (p.dark_green(i, :) + p.black(i, :)) / 2;
end

% Interpolate between the two color schemes
p.dark_blue = zeros(num_colors, 3);

for i = 1:num_colors
    p.dark_blue(i, :) = (p.blue(i, :) + p.black(i, :)) / 2;
end

% Interpolate between the two color schemes
p.darker_blue = zeros(num_colors, 3);

for i = 1:num_colors
    p.darker_blue(i, :) = (p.dark_blue(i, :) + p.black(i, :)) / 2;
end



end