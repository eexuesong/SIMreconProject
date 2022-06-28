clc;
clear all;
close all;

x_offset = 0;   % Screen position
y_offset = 0;   % Screen position
width  = 512; % Width of figure
height = 512; % Height of figure (by default in pixels)
fill_in_SIM_gap = true;

%% Physical parameters
NA = 1.35;
pixel_size = 40;     % 40 nm
k_max = 1 / (2 * pixel_size);
lambda_Ex = 488;    % unit: nanometer
lambda_Em = 525;
RI = 1.406;
pupil_filling_factor = 0.92;

%% Emission OTF
krmax = 2 * NA / lambda_Em;         % resolution =  lambda_Em / (2 * NA); krmax = 1 / resolution;
krmax_half = NA / lambda_Em;
radius = RI / lambda_Em;                 % radius of wide-field OTF torus;
radius_I2M = 2 * radius;                       % radius of I2M OTF torus;

%% Excitation OTF
kr = NA * pupil_filling_factor / lambda_Ex;
Beta = asin(NA * pupil_filling_factor / RI);
kz1 = (1 - cos(Beta)) * RI / lambda_Ex;
kz2 = 2 * kz1;
kz3 = (1 + cos(Beta)) * RI / lambda_Ex;
kz4 = 2 * RI / lambda_Ex;

%% Scatter Dots coordinates
WF_r = 0;
WF_z = 0;
SW_r = [0, 0];
SW_z = [kz4, -kz4];
SIM_r = [-2 * kr, -kr, -kr, 0, kr, kr, 2 * kr];
SIM_z = [0, kz1, -kz1, 0, kz1, -kz1, 0];
SW_SIM_r = [-kr, -kr, 0, 0, kr, kr];
SW_SIM_z = [kz3, -kz3, kz4, -kz4, kz3, -kz3];
I5S_r = [-2 * kr, -2 * kr, -kr, -kr, 0, 0, 0, 0, kr, kr, 2 * kr, 2 * kr];
I5S_z = [kz2, -kz2, kz3, -kz3, kz4, kz2, -kz2, -kz4, kz3, -kz3, kz2, -kz2];

%% Generate 2D wide-field OTF profiles
% phi is the counterclockwise angle in the r-z plane measured in radians from the positive z-axis
phi_max = asin(NA / RI);    % phi_max is the maximum phi angle corresponding to the emmision OTF when kr = 0
phi_min = -phi_max;

r = (-krmax_half : krmax_half/100 : krmax_half);
z_offset = radius * cos(phi_max);
z_top = sqrt(radius ^ 2 - r.^2) - z_offset;
z_bottom = -z_top;
% z_highest = radius - z_offset = radius * (1 - cos(phi_max));
z_peak = radius * (1 - cos(phi_max));

% order = 0
r_left = r - krmax_half;
r_right = r + krmax_half;

% order = 1, positive direction
r_left_positive = r_left + kr;
r_right_positive = r_right + kr;
% order = 1, negative direction
r_left_negative = r_left - kr;
r_right_negative = r_right - kr;

% order = 2, positive direction
r_left_positive_2 = r_left + 2 * kr;
r_right_positive_2 = r_right + 2 *kr;
% order = 2, negative direction
r_left_negative_2 = r_left - 2 * kr;
r_right_negative_2 = r_right - 2 * kr;

%% Generate 2D I2M OTF profiles
% phi is the counterclockwise angle in the r-z plane measured in radians from the positive z-axis
phi_I2M_max = asin(krmax / radius_I2M);    % phi_I2M_max is the maximum phi angle corresponding to the emmision OTF when kr = 0
phi_I2M_min = -phi_I2M_max;

% order = 0;
r_I2M = (-krmax : krmax/200 : krmax);
z_I2M_top = sqrt(radius_I2M ^ 2 - r_I2M.^2);
z_I2M_line_top = radius_I2M * cos(phi_I2M_max) * ones(size(r_I2M));
z_I2M_bottom = - sqrt(radius_I2M ^ 2 - r_I2M.^2);
z_I2M_line_bottom = -radius_I2M * cos(phi_I2M_max) * ones(size(r_I2M));

% order = 1, positive direction
r_I2M_positive = r_I2M + kr;
% order = 1, negative direction
r_I2M_negative = r_I2M - kr;

% order = 2, positive direction
r_I2M_positive_2 = r_I2M + 2 * kr;
% order = 2, negative direction
r_I2M_negative_2 = r_I2M - 2 * kr;


%% Wide-field OTF
fig = figure('Name', 'WF', 'Position', [x_offset y_offset width height], 'Units', 'pixels');
hold all;
% patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% patch([r_right, fliplr(r_right)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
scatter(WF_r, WF_z, 'filled','r')

% arrow_horizontal_x = [0.098, 0.2];   % adjust length and location of arrow 
% arrow_horizontal_y = [0.1, 0.1];
% annotation('textarrow', arrow_horizontal_x, arrow_horizontal_y, 'FontSize', 13, 'Linewidth',2);
% dim = [0.12, 0.05, 0.05, 0.05];  %  Size and location: four-element vector of the form [x y w h]
% t = annotation('textbox', dim);
% t.String = 'k_{r}';
% t.FontSize = 13;
% t.FontAngle = 'italic';
% t.FontWeight = 'bold';
% t.EdgeColor = 'none';
% t.LineWidth = 2;
% t.VerticalAlignment = 'top';

% arrow_vertical_x = [0.1, 0.1];  % adjust length and location of arrow 
% arrow_vertical_y = [0.1, 0.2];  
% annotation('textarrow', arrow_vertical_x, arrow_vertical_y, 'FontSize', 13, 'Linewidth',2);
% arrow_text_x = [0.065, 0.165];
% arrow_text_y = [0.165, 0.165];
% t = annotation('textarrow', arrow_text_x, arrow_text_y);
% t.String = 'k_{z}';
% t.FontSize = 13;
% t.FontAngle = 'italic';
% t.FontWeight = 'bold';
% t.Color = 'k';
% t.HeadStyle = 'none';
% t.LineStyle = 'none';
% t.TextRotation = 90;
% 
% scale_bar_dim = [0.8, 0.1, 0.1, 0]; % Size and location, specified as a four-element vector of the form [x_begin y_begin dx dy].
% h = annotation('line');
% h.Position = scale_bar_dim;
% h.Color = 'k';
% h.LineStyle = '-';
% h.LineWidth = 2;
% dim = [0.77, 0.1, 0.2, 0.05];  %  Size and location: four-element vector of the form [x y w h]
% t = annotation('textbox', dim);
% t.String = strcat('^{1}/_{', num2str(10 * pixel_size / 2),'} nm^{-1}');
% t.FontName = 'Serif';
% t.FontSize = 11;
% % t.FontAngle = 'italic';
% t.FontWeight = 'bold';
% t.EdgeColor = 'none';
% t.LineWidth = 2;
% t.VerticalAlignment = 'bottom';

hold off;
axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);
% saveas(fig, 'WF', 'tiffn');

%% I2M OTF
figure('Name', 'I2M', 'Position', [x_offset y_offset width height]);
hold all;
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_I2M, fliplr(r_I2M)], [z_I2M_top, fliplr(z_I2M_line_top)], 'b', 'EdgeColor', 'none');
patch([r_I2M, fliplr(r_I2M)], [z_I2M_bottom, fliplr(z_I2M_line_bottom)], 'b', 'EdgeColor', 'none');
scatter(WF_r, WF_z, 'filled','r')
hold off;

axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);


%% SW OTF
figure('Name', 'SW', 'Position', [x_offset y_offset width height]);
hold all;

% WF
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
% SW: kz4
for index = 0:1
    z_top_4 = z_top + (-1)^index * kz4;
    z_bottom_4 = z_bottom + (-1)^index * kz4;
    patch([r_left, fliplr(r_left)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
    patch([r_right, fliplr(r_right)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
end

scatter(WF_r, WF_z, 'filled','r')
hold on;
scatter(SW_r, SW_z, 'filled','r')
hold off;
axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);


%% 3D-SIM OTF
figure('Position', [x_offset y_offset width height]);
hold all;
%% order = 0
% WF
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_bottom, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');

%% order = 1, positive direction;
% 3D-SIM: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_positive, fliplr(r_left_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive, fliplr(r_right_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
end
%% order = 1, negative direction;
% 3D-SIM: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_negative, fliplr(r_left_negative)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative, fliplr(r_right_negative)], [z_bottom_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
end
if fill_in_SIM_gap
    for index = 0:1
        r_filled_in = r + (-1)^index * kr;
        z_top_filled_in = (z_peak + kz1) * ones(size(r_filled_in));
        z_bottom_filled_in = kz1 * ones(size(r_filled_in));        
        patch([r_filled_in, fliplr(r_filled_in)], [z_top_filled_in, fliplr(z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
        patch([r_filled_in, fliplr(r_filled_in)], [-z_top_filled_in, fliplr(-z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
    end
end
%% order = 2, positive direction;
patch([r_left_positive_2, fliplr(r_left_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_positive_2, fliplr(r_right_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
%% order = 2, negative direction;
patch([r_left_negative_2, fliplr(r_left_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_negative_2, fliplr(r_right_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');

scatter(SIM_r, SIM_z, 'filled','r')
hold off;
axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);


%% SW-SIM OTF
figure('Name', 'SW SIM', 'Position', [x_offset y_offset width height]);
hold all;
%% order = 0
% WF
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_bottom, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
% SW: kz4
for index = 0:1
    z_top_4 = z_top + (-1)^index * kz4;
    z_bottom_4 = z_bottom + (-1)^index * kz4;
    patch([r_left, fliplr(r_left)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
    patch([r_right, fliplr(r_right)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
end
%% order = 1, positive direction;
% 3D-SIM: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_positive, fliplr(r_left_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive, fliplr(r_right_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
end
% 3D-SIM: kz3
for index = 0:1
    z_top_3 = z_top + (-1)^index * kz3;
    z_bottom_3 = z_bottom + (-1)^index * kz3;
    patch([r_left_positive, fliplr(r_left_positive)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive, fliplr(r_right_positive)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
end
%% order = 1, negative direction;
% 3D-SIM: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_negative, fliplr(r_left_negative)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative, fliplr(r_right_negative)], [z_bottom_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
end
% 3D-SIM: kz3
for index = 0:1
    z_top_3 = z_top + (-1)^index * kz3;
    z_bottom_3 = z_bottom + (-1)^index * kz3;
    patch([r_left_negative, fliplr(r_left_negative)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative, fliplr(r_right_negative)], [z_bottom_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
end
if fill_in_SIM_gap
    for index = 0:1
        r_filled_in = r + (-1)^index * kr;
        z_top_filled_in = (z_peak + kz1) * ones(size(r_filled_in));
        z_bottom_filled_in = kz1 * ones(size(r_filled_in));        
        patch([r_filled_in, fliplr(r_filled_in)], [z_top_filled_in, fliplr(z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
        patch([r_filled_in, fliplr(r_filled_in)], [-z_top_filled_in, fliplr(-z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
    end
    
    for index = 0:1
        r_filled_in = r + (-1)^index * kr;
        z_top_filled_in = kz3 * ones(size(r_filled_in));  
        z_bottom_filled_in = (kz3 - z_peak) * ones(size(r_filled_in));  
        patch([r_filled_in, fliplr(r_filled_in)], [z_top_filled_in, fliplr(z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
        patch([r_filled_in, fliplr(r_filled_in)], [-z_top_filled_in, fliplr(-z_bottom_filled_in)], 'b', 'EdgeColor', 'none');
    end
end
%% order = 2, positive direction;
patch([r_left_positive_2, fliplr(r_left_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_positive_2, fliplr(r_right_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
%% order = 2, negative direction;
patch([r_left_negative_2, fliplr(r_left_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_negative_2, fliplr(r_right_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');

scatter(SIM_r, SIM_z, 'filled','r')
scatter(SW_SIM_r, SW_SIM_z, 'filled','g')

%% Scale bar
arrow_horizontal_x = [0.098, 0.2];   % adjust length and location of arrow 
arrow_horizontal_y = [0.1, 0.1];
annotation('textarrow', arrow_horizontal_x, arrow_horizontal_y, 'FontSize', 13, 'Linewidth',2);
dim = [0.12, 0.05, 0.05, 0.05];  %  Size and location: four-element vector of the form [x y w h]
t = annotation('textbox', dim);
t.String = 'k_{r}';
t.FontSize = 20;
t.FontAngle = 'italic';
t.FontWeight = 'bold';
t.EdgeColor = 'none';
t.LineWidth = 2;
t.VerticalAlignment = 'top';

arrow_vertical_x = [0.1, 0.1];  % adjust length and location of arrow 
arrow_vertical_y = [0.1, 0.2];  
annotation('textarrow', arrow_vertical_x, arrow_vertical_y, 'FontSize', 13, 'Linewidth',2);
arrow_text_x = [0.05, 0.165];
arrow_text_y = [0.165, 0.165];
t = annotation('textarrow', arrow_text_x, arrow_text_y);
t.String = 'k_{z}';
t.FontSize = 20;
t.FontAngle = 'italic';
t.FontWeight = 'bold';
t.Color = 'k';
t.HeadStyle = 'none';
t.LineStyle = 'none';
t.TextRotation = 90;

scale_bar_dim = [0.8, 0.1, 0.1, 0]; % Size and location, specified as a four-element vector of the form [x_begin y_begin dx dy].
h = annotation('line');
h.Position = scale_bar_dim;
h.Color = 'k';
h.LineStyle = '-';
h.LineWidth = 2;
dim = [0.75, 0.1, 0.4, 0.1];  %  Size and location: four-element vector of the form [x y w h]
t = annotation('textbox', dim);
t.String = strcat('^{1}/_{', num2str(10 * pixel_size / 2),'} nm^{-1}');
t.FontName = 'Serif';
t.FontSize = 20;
% t.FontAngle = 'italic';
t.FontWeight = 'bold';
t.EdgeColor = 'none';
t.LineWidth = 2;
t.VerticalAlignment = 'bottom';

hold off;
axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);


%% I5S OTF
figure('Name', 'I5S', 'Position', [x_offset y_offset width height]);
hold all;
%% order = 0
% WF
patch([r_left, fliplr(r_left)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right, fliplr(r_right)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_I2M, fliplr(r_I2M)], [z_I2M_top, fliplr(z_I2M_line_top)], 'b', 'EdgeColor', 'none');
patch([r_I2M, fliplr(r_I2M)], [z_I2M_bottom, fliplr(z_I2M_line_bottom)], 'b', 'EdgeColor', 'none');
% SW: kz4
for index = 0:1
    z_top_4 = z_top + (-1)^index * kz4;
    z_bottom_4 = z_bottom + (-1)^index * kz4;
    patch([r_left, fliplr(r_left)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
    patch([r_right, fliplr(r_right)], [z_top_4, fliplr(z_bottom_4)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_4 = z_I2M_top + (-1)^index * kz4;
    z_I2M_line_top_4 = z_I2M_line_top + (-1)^index * kz4;
    z_I2M_bottom_4 = z_I2M_bottom + (-1)^index * kz4;
    z_I2M_line_bottom_4 = z_I2M_line_bottom + (-1)^index * kz4;
    patch([r_I2M, fliplr(r_I2M)], [z_I2M_top_4, fliplr(z_I2M_line_top_4)], 'b', 'EdgeColor', 'none');
    patch([r_I2M, fliplr(r_I2M)], [z_I2M_bottom_4, fliplr(z_I2M_line_bottom_4)], 'b', 'EdgeColor', 'none');
end

% SW: kz2
for index = 0:1
    z_top_2 = z_top + (-1)^index * kz2;
    z_bottom_2 = z_bottom + (-1)^index * kz2;
    patch([r_left, fliplr(r_left)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    patch([r_right, fliplr(r_right)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_2 = z_I2M_top + (-1)^index * kz2;
    z_I2M_line_top_2 = z_I2M_line_top + (-1)^index * kz2;
    z_I2M_bottom_2 = z_I2M_bottom + (-1)^index * kz2;
    z_I2M_line_bottom_2 = z_I2M_line_bottom + (-1)^index * kz2;
    patch([r_I2M, fliplr(r_I2M)], [z_I2M_top_2, fliplr(z_I2M_line_top_2)], 'b', 'EdgeColor', 'none');
    patch([r_I2M, fliplr(r_I2M)], [z_I2M_bottom_2, fliplr(z_I2M_line_bottom_2)], 'b', 'EdgeColor', 'none');
end
%% order = 1, positive direction;
% SW: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_positive, fliplr(r_left_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive, fliplr(r_right_positive)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_1 = z_I2M_top + (-1)^index * kz1;
    z_I2M_line_top_1 = z_I2M_line_top + (-1)^index * kz1;
    z_I2M_bottom_1 = z_I2M_bottom + (-1)^index * kz1;
    z_I2M_line_bottom_1 = z_I2M_line_bottom + (-1)^index * kz1;
    patch([r_I2M_positive, fliplr(r_I2M_positive)], [z_I2M_top_1, fliplr(z_I2M_line_top_1)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_positive, fliplr(r_I2M_positive)], [z_I2M_bottom_1, fliplr(z_I2M_line_bottom_1)], 'b', 'EdgeColor', 'none');
end

% SW: kz3
for index = 0:1
    z_top_3 = z_top + (-1)^index * kz3;
    z_bottom_3 = z_bottom + (-1)^index * kz3;
    patch([r_left_positive, fliplr(r_left_positive)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive, fliplr(r_right_positive)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_3 = z_I2M_top + (-1)^index * kz3;
    z_I2M_line_top_3 = z_I2M_line_top + (-1)^index * kz3;
    z_I2M_bottom_3 = z_I2M_bottom + (-1)^index * kz3;
    z_I2M_line_bottom_3 = z_I2M_line_bottom + (-1)^index * kz3;
    patch([r_I2M_positive, fliplr(r_I2M_positive)], [z_I2M_top_3, fliplr(z_I2M_line_top_3)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_positive, fliplr(r_I2M_positive)], [z_I2M_bottom_3, fliplr(z_I2M_line_bottom_3)], 'b', 'EdgeColor', 'none');
end
%% order = 1, negative direction;
% SW: kz1
for index = 0:1
    z_top_1 = z_top + (-1)^index * kz1;
    z_bottom_1 = z_bottom + (-1)^index * kz1;
    patch([r_left_negative, fliplr(r_left_negative)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative, fliplr(r_right_negative)], [z_top_1, fliplr(z_bottom_1)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_1 = z_I2M_top + (-1)^index * kz1;
    z_I2M_line_top_1 = z_I2M_line_top + (-1)^index * kz1;
    z_I2M_bottom_1 = z_I2M_bottom + (-1)^index * kz1;
    z_I2M_line_bottom_1 = z_I2M_line_bottom + (-1)^index * kz1;
    patch([r_I2M_negative, fliplr(r_I2M_negative)], [z_I2M_top_1, fliplr(z_I2M_line_top_1)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_negative, fliplr(r_I2M_negative)], [z_I2M_bottom_1, fliplr(z_I2M_line_bottom_1)], 'b', 'EdgeColor', 'none');
end

% 3D-SIM: kz3
for index = 0:1
    z_top_3 = z_top + (-1)^index * kz3;
    z_bottom_3 = z_bottom + (-1)^index * kz3;
    patch([r_left_negative, fliplr(r_left_negative)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative, fliplr(r_right_negative)], [z_top_3, fliplr(z_bottom_3)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_3 = z_I2M_top + (-1)^index * kz3;
    z_I2M_line_top_3 = z_I2M_line_top + (-1)^index * kz3;
    z_I2M_bottom_3 = z_I2M_bottom + (-1)^index * kz3;
    z_I2M_line_bottom_3 = z_I2M_line_bottom + (-1)^index * kz3;
    patch([r_I2M_negative, fliplr(r_I2M_negative)], [z_I2M_top_3, fliplr(z_I2M_line_top_3)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_negative, fliplr(r_I2M_negative)], [z_I2M_bottom_3, fliplr(z_I2M_line_bottom_3)], 'b', 'EdgeColor', 'none');
end
%% order = 2, positive direction;
patch([r_left_positive_2, fliplr(r_left_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_positive_2, fliplr(r_right_positive_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_I2M_positive_2, fliplr(r_I2M_positive_2)], [z_I2M_top, fliplr(z_I2M_line_top)], 'b', 'EdgeColor', 'none');
patch([r_I2M_positive_2, fliplr(r_I2M_positive_2)], [z_I2M_bottom, fliplr(z_I2M_line_bottom)], 'b', 'EdgeColor', 'none');

% SW: kz2
for index = 0:1
    z_top_2 = z_top + (-1)^index * kz2;
    z_bottom_2 = z_bottom + (-1)^index * kz2;
    patch([r_left_positive_2, fliplr(r_left_positive_2)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    patch([r_right_positive_2, fliplr(r_right_positive_2)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_2 = z_I2M_top + (-1)^index * kz2;
    z_I2M_line_top_2 = z_I2M_line_top + (-1)^index * kz2;
    z_I2M_bottom_2 = z_I2M_bottom + (-1)^index * kz2;
    z_I2M_line_bottom_2 = z_I2M_line_bottom + (-1)^index * kz2;
    patch([r_I2M_positive_2, fliplr(r_I2M_positive_2)], [z_I2M_top_2, fliplr(z_I2M_line_top_2)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_positive_2, fliplr(r_I2M_positive_2)], [z_I2M_bottom_2, fliplr(z_I2M_line_bottom_2)], 'b', 'EdgeColor', 'none');
end
%% order = 2, negative direction;
patch([r_left_negative_2, fliplr(r_left_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_right_negative_2, fliplr(r_right_negative_2)], [z_top, fliplr(z_bottom)], 'b', 'EdgeColor', 'none');
patch([r_I2M_negative_2, fliplr(r_I2M_negative_2)], [z_I2M_top, fliplr(z_I2M_line_top)], 'b', 'EdgeColor', 'none');
patch([r_I2M_negative_2, fliplr(r_I2M_negative_2)], [z_I2M_bottom, fliplr(z_I2M_line_bottom)], 'b', 'EdgeColor', 'none');

% SW: kz2
for index = 0:1
    z_top_2 = z_top + (-1)^index * kz2;
    z_bottom_2 = z_bottom + (-1)^index * kz2;
    patch([r_left_negative_2, fliplr(r_left_negative_2)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    patch([r_right_negative_2, fliplr(r_right_negative_2)], [z_top_2, fliplr(z_bottom_2)], 'b', 'EdgeColor', 'none');
    
    z_I2M_top_2 = z_I2M_top + (-1)^index * kz2;
    z_I2M_line_top_2 = z_I2M_line_top + (-1)^index * kz2;
    z_I2M_bottom_2 = z_I2M_bottom + (-1)^index * kz2;
    z_I2M_line_bottom_2 = z_I2M_line_bottom + (-1)^index * kz2;
    patch([r_I2M_negative_2, fliplr(r_I2M_negative_2)], [z_I2M_top_2, fliplr(z_I2M_line_top_2)], 'b', 'EdgeColor', 'none');
    patch([r_I2M_negative_2, fliplr(r_I2M_negative_2)], [z_I2M_bottom_2, fliplr(z_I2M_line_bottom_2)], 'b', 'EdgeColor', 'none');
end

scatter(SIM_r, SIM_z, 'filled','r')
scatter(I5S_r, I5S_z, 'filled','g')

hold off;
axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max]);
