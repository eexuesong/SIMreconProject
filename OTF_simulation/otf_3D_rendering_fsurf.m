clc;
clear all;
close all;

x_offset = 0;   % Screen position
y_offset = 0;   % Screen position
width  = 512; % Width of figure
height = 512; % Height of figure (by default in pixels)

NA = 1.35;
pixel_size = 40;    % 40 nm
k_max = 1 / (2 * pixel_size);
lambda_Ex = 488;    % unit: nanometer
lambda_Em = 525;
RI = 1.406;
pupil_filling_factor = 0.92;
angle = [0, pi/3, 2*pi/3];

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


%% Generate 3D-OTF profiles
% theta is the counterclockwise angle in the x-y plane measured in radians from the positive x-axis
% phi is the counterclockwise angle in the x-z plane measured in radians from the positive z-axis
phi_max = asin(NA / RI);    % phi_max is the maximum phi angle corresponding to the emmision OTF when kr = 0
phi_min = -phi_max;
theta_min = 0;
% theta_max = pi;             % Cross-sectional view
theta_max = 2 * pi;          % Full view

% order = 0
% kr = krmax_half - r * sin(phi) = NA / lambda_Em - RI / lambda_Em * sin(phi) = (NA - RI * sin(phi) / lambda_Em
% x = kr * cos(theta);
% y = kr * sin(theta);
x = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta);
y = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta);
z_top = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max));
z_bottom = @(theta, phi) -(RI / lambda_Em) * (cos(phi) - cos(phi_max));


%% Generate 3D I2M OTF profiles
% phi_I2M is the counterclockwise angle in the r-z plane measured in radians from the positive z-axis
phi_I2M_max = asin(krmax / radius_I2M);    % phi_I2M_max is the maximum phi angle corresponding to the emmision OTF when kr = 0
phi_I2M_min = -phi_I2M_max;
theta_I2M_min = 0;
% theta_max = pi;             % Cross-sectional view
theta_I2M_max = 2 * pi;          % Full view

% order = 0;
x_I2M = @(theta_I2M, phi_I2M) radius_I2M * sin(phi_I2M) * cos(theta_I2M);
y_I2M = @(theta_I2M, phi_I2M) radius_I2M * sin(phi_I2M) * sin(theta_I2M);
z_I2M_top = @(theta_I2M, phi_I2M) radius_I2M * cos(phi_I2M);
z_I2M_bottom = @(theta_I2M, phi_I2M) -radius_I2M * cos(phi_I2M);
z_I2M_line_top = @(theta_I2M, phi_I2M)radius_I2M * cos(phi_I2M_max);
z_I2M_line_bottom = @(theta_I2M, phi_I2M) -radius_I2M * cos(phi_I2M_max);


% %% Wide-field OTF
% figure('Name', 'WF', 'Position', [x_offset y_offset width height]);
% 
% h_top = fsurf(x, y, z_top, [theta_min, theta_max, phi_min, phi_max]);
% % h_top.FaceColor = '#808080';
% h_top.FaceColor = [0 0.4470 0.7410];
% h_top.FaceAlpha = 1;    % A value of 1 is fully opaque and 0 is completely transparent. Values between 0 and 1 are semitransparent.
% h_top.EdgeColor = 'none';
% h_top.AmbientStrength = 0.3;
% % h_top.SpecularColorReflectance = 0.5;
% hold all;
% 
% h_bottom = fsurf(x, y, z_bottom, [theta_min, theta_max, phi_min, phi_max]);
% h_bottom.FaceColor = [0 0.4470 0.7410];
% h_bottom.FaceAlpha = 1;
% h_bottom.EdgeColor = 'none';
% h_bottom.AmbientStrength = 0.3;
% hold off;
% 
% axis equal;
% axis off;
% axis([-k_max, k_max, -k_max, k_max, -k_max, k_max]);
% camlight('right');
% 
% %% I2M OTF
% figure('Name', 'I2M', 'Position', [x_offset y_offset width height]);
% 
% h_top = fsurf(x, y, z_top, [theta_min, theta_max, phi_min, phi_max]);
% % h_top.FaceColor = '#808080';
% h_top.FaceColor = [0 0.4470 0.7410];
% h_top.FaceAlpha = 1;    % A value of 1 is fully opaque and 0 is completely transparent. Values between 0 and 1 are semitransparent.
% h_top.EdgeColor = 'none';
% h_top.AmbientStrength = 0.3;
% % h_top.SpecularColorReflectance = 0.5;
% hold all;
% 
% h_bottom = fsurf(x, y, z_bottom, [theta_min, theta_max, phi_min, phi_max]);
% h_bottom.FaceColor = [0 0.4470 0.7410];
% h_bottom.FaceAlpha = 1;
% h_bottom.EdgeColor = 'none';
% h_bottom.AmbientStrength = 0.3;
% 
% h_top_I2M = fsurf(x_I2M, y_I2M, z_I2M_top, [theta_I2M_min, theta_I2M_max, phi_I2M_min, phi_I2M_max]);
% h_top_I2M.FaceColor = [0 0.4470 0.7410];
% h_top_I2M.FaceAlpha = 1;
% h_top_I2M.EdgeColor = 'none';
% h_top_I2M.AmbientStrength = 0.3;
% 
% h_top_I2M = fsurf(x_I2M, y_I2M, z_I2M_line_top, [theta_I2M_min, theta_I2M_max, phi_I2M_min, phi_I2M_max]);
% h_top_I2M.FaceColor = [0 0.4470 0.7410];
% h_top_I2M.FaceAlpha = 1;
% h_top_I2M.EdgeColor = 'none';
% h_top_I2M.AmbientStrength = 0.3;
% 
% h_bottom_I2M = fsurf(x_I2M, y_I2M, z_I2M_bottom, [theta_I2M_min, theta_I2M_max, phi_I2M_min, phi_I2M_max]);
% h_bottom_I2M.FaceColor = [0 0.4470 0.7410];
% h_bottom_I2M.FaceAlpha = 1;
% h_bottom_I2M.EdgeColor = 'none';
% h_bottom_I2M.AmbientStrength = 0.3;
% 
% h_bottom_I2M = fsurf(x_I2M, y_I2M, z_I2M_line_bottom, [theta_I2M_min, theta_I2M_max, phi_I2M_min, phi_I2M_max]);
% h_bottom_I2M.FaceColor = [0 0.4470 0.7410];
% h_bottom_I2M.FaceAlpha = 1;
% h_bottom_I2M.EdgeColor = 'none';
% h_bottom_I2M.AmbientStrength = 0.3;
% 
% axis equal;
% axis off;
% axis([-k_max, k_max, -k_max, k_max, -k_max, k_max]);
% camlight('right');
% 
% 
% %% SW OTF
% figure('Name', 'SW', 'Position', [x_offset y_offset width height]);
% 
% % WF
% h_top = fsurf(x, y, z_top, [theta_min, theta_max, phi_min, phi_max]);
% h_top.FaceColor = [0 0.4470 0.7410];
% h_top.FaceAlpha = 1;
% h_top.EdgeColor = 'none';
% h_top.AmbientStrength = 0.3;
% hold all;
% 
% h_bottom = fsurf(x, y, z_bottom, [theta_min, theta_max, phi_min, phi_max]);
% h_bottom.FaceColor = [0 0.4470 0.7410];
% h_bottom.FaceAlpha = 1;
% h_bottom.EdgeColor = 'none';
% h_bottom.AmbientStrength = 0.3;
% 
% % SW: kz4
% for index = 0:1
%     z_top_4 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz4;
%     z_bottom_4 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz4;
%     h_top = fsurf(x, y, z_top_4, [theta_min, theta_max, phi_min, phi_max]);
%     h_bottom = fsurf(x, y, z_bottom_4, [theta_min, theta_max, phi_min, phi_max]);
%     
%     h_top.FaceColor = [0 0.4470 0.7410];
%     h_top.FaceAlpha = 1;
%     h_top.EdgeColor = 'none';
%     h_top.AmbientStrength = 0.3;
%     
%     h_bottom.FaceColor = [0 0.4470 0.7410];
%     h_bottom.FaceAlpha = 1;
%     h_bottom.EdgeColor = 'none';
%     h_bottom.AmbientStrength = 0.3;
% end
% hold off;
% 
% axis equal;
% axis off;
% axis([-k_max, k_max, -k_max, k_max, -k_max, k_max]);
% camlight('right');
% 
% %% 3D-SIM OTF (3 directions)
% figure('Name', '3D SIM', 'Position', [x_offset y_offset width height]);
% %% order = 0
% % WF
% h_top = fsurf(x, y, z_top, [theta_min, theta_max, phi_min, phi_max]);
% h_top.FaceColor = [0 0.4470 0.7410];
% h_top.FaceAlpha = 1;
% h_top.EdgeColor = 'none';
% h_top.AmbientStrength = 0.3;
% hold all;
% 
% h_bottom = fsurf(x, y, z_bottom, [theta_min, theta_max, phi_min, phi_max]);
% h_bottom.FaceColor = [0 0.4470 0.7410];
% h_bottom.FaceAlpha = 1;
% h_bottom.EdgeColor = 'none';
% h_bottom.AmbientStrength = 0.3;
% 
% for dir = 1:3
% 
%     %% order = 1, positive direction;
%     % 3D-SIM: kz1
%     x_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) + kr * cos(angle(dir));
%     y_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) + kr * sin(angle(dir));
%     
%     for index = 0:1
%         z_top_1 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
%         z_bottom_1 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
%         h_top = fsurf(x_positive, y_positive, z_top_1, [theta_min, theta_max, phi_min, phi_max]);
%         h_bottom = fsurf(x_positive, y_positive, z_bottom_1, [theta_min, theta_max, phi_min, phi_max]);
%         
%         h_top.FaceColor = [0 0.4470 0.7410];
%         h_top.FaceAlpha = 1;
%         h_top.EdgeColor = 'none';
%         h_top.AmbientStrength = 0.3;
%         
%         h_bottom.FaceColor = [0 0.4470 0.7410];
%         h_bottom.FaceAlpha = 1;
%         h_bottom.EdgeColor = 'none';
%         h_bottom.AmbientStrength = 0.3;
%     end
%     
%     %% order = 1, nagetive direction;
%     % 3D-SIM: kz1
%     x_negative = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) - kr * cos(angle(dir));
%     y_negative = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) - kr * sin(angle(dir));
%     
%     for index = 0:1
%         z_top_1 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
%         z_bottom_1 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
%         h_top = fsurf(x_negative, y_negative, z_top_1, [theta_min, theta_max, phi_min, phi_max]);
%         h_bottom = fsurf(x_negative, y_negative, z_bottom_1, [theta_min, theta_max, phi_min, phi_max]);
%         
%         h_top.FaceColor = [0 0.4470 0.7410];
%         h_top.FaceAlpha = 1;
%         h_top.EdgeColor = 'none';
%         h_top.AmbientStrength = 0.3;
%     
%         h_bottom.FaceColor = [0 0.4470 0.7410];
%         h_bottom.FaceAlpha = 1;
%         h_bottom.EdgeColor = 'none';
%         h_bottom.AmbientStrength = 0.3;
%     end
%     
%     %% order = 2
%     z_top = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max));
%     z_bottom = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max));
%     % positive direction
%     x_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) + 2 * kr * cos(angle(dir));
%     y_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) + 2 * kr * sin(angle(dir));
%     h_top = fsurf(x_positive, y_positive, z_top,[theta_min, theta_max, phi_min, phi_max]);
%     h_bottom = fsurf(x_positive, y_positive, z_bottom,[theta_min, theta_max, phi_min, phi_max]);
%     
%     h_top.FaceColor = [0 0.4470 0.7410];
%     h_top.FaceAlpha = 1;
%     h_top.EdgeColor = 'none';
%     h_top.AmbientStrength = 0.3;
% 
%     h_bottom.FaceColor = [0 0.4470 0.7410];
%     h_bottom.FaceAlpha = 1;
%     h_bottom.EdgeColor = 'none';
%     h_bottom.AmbientStrength = 0.3;
% 
%     % negative direction
%     x_negetive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) - 2 * kr * cos(angle(dir));
%     y_negetive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) - 2 * kr * sin(angle(dir));
%     h_top = fsurf(x_negetive, y_negetive, z_top,[theta_min, theta_max, phi_min, phi_max]);
%     h_bottom = fsurf(x_negetive, y_negetive, z_bottom,[theta_min, theta_max, phi_min, phi_max]);
%     
%     h_top.FaceColor = [0 0.4470 0.7410];
%     h_top.FaceAlpha = 1;
%     h_top.EdgeColor = 'none';
%     h_top.AmbientStrength = 0.3;
%     
%     h_bottom.FaceColor = [0 0.4470 0.7410];
%     h_bottom.FaceAlpha = 1;
%     h_bottom.EdgeColor = 'none';
%     h_bottom.AmbientStrength = 0.3;
% end
% hold off
% 
% axis equal;
% axis off;
% axis([-k_max, k_max, -k_max, k_max, -k_max, k_max]);
% camlight('right');


%% SW-SIM OTF (3 directions)
figure('Name', 'SW SIM', 'Position', [x_offset y_offset width height]);
%% order = 0
% WF
h_top = fsurf(x, y, z_top,[theta_min, theta_max, phi_min, phi_max]);
h_top.FaceColor = [0 0.4470 0.7410];
h_top.FaceAlpha = 1;
h_top.EdgeColor = 'none';
h_top.AmbientStrength = 0.3;
hold all;

h_bottom = fsurf(x, y, z_bottom,[theta_min, theta_max, phi_min, phi_max]);
h_bottom.FaceColor = [0 0.4470 0.7410];
h_bottom.FaceAlpha = 1;
h_bottom.EdgeColor = 'none';
h_bottom.AmbientStrength = 0.3;

% SW: kz4
for index = 0:1
    z_top_4 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz4;
    z_bottom_4 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz4;
    h_top = fsurf(x, y, z_top_4, [theta_min, theta_max, phi_min, phi_max]);
    h_bottom = fsurf(x, y, z_bottom_4, [theta_min, theta_max, phi_min, phi_max]);
    
    h_top.FaceColor = [0 0.4470 0.7410];
    h_top.FaceAlpha = 1;
    h_top.EdgeColor = 'none';
    h_top.AmbientStrength = 0.3;
    
    h_bottom.FaceColor = [0 0.4470 0.7410];
    h_bottom.FaceAlpha = 1;
    h_bottom.EdgeColor = 'none';
    h_bottom.AmbientStrength = 0.3;
end

for dir = 1:3
    %% order = 1, positive direction;
    x_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) + kr * cos(angle(dir));
    y_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) + kr * sin(angle(dir));
    
    % 3D-SIM: kz1
    for index = 0:1
        z_top_1 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
        z_bottom_1 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
        h_top = fsurf(x_positive, y_positive, z_top_1, [theta_min, theta_max, phi_min, phi_max]);
        h_bottom = fsurf(x_positive, y_positive, z_bottom_1, [theta_min, theta_max, phi_min, phi_max]);
        
        h_top.FaceColor = [0 0.4470 0.7410];
        h_top.FaceAlpha = 1;
        h_top.EdgeColor = 'none';
        h_top.AmbientStrength = 0.3;
        
        h_bottom.FaceColor = [0 0.4470 0.7410];
        h_bottom.FaceAlpha = 1;
        h_bottom.EdgeColor = 'none';
        h_bottom.AmbientStrength = 0.3;
    end
    
    % 3D-SIM: kz3
    for index = 0:1
        z_top_3 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz3;
        z_bottom_3 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz3;
        h_top = fsurf(x_positive, y_positive, z_top_3, [theta_min, theta_max, phi_min, phi_max]);
        h_bottom = fsurf(x_positive, y_positive, z_bottom_3, [theta_min, theta_max, phi_min, phi_max]);
        
        h_top.FaceColor = [0 0.4470 0.7410];
        h_top.FaceAlpha = 1;
        h_top.EdgeColor = 'none';
        h_top.AmbientStrength = 0.3;
        
        h_bottom.FaceColor = [0 0.4470 0.7410];
        h_bottom.FaceAlpha = 1;
        h_bottom.EdgeColor = 'none';
        h_bottom.AmbientStrength = 0.3;
    end
    
    %% order = 1, nagetive direction;
    x_negative = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) - kr * cos(angle(dir));
    y_negative = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) - kr * sin(angle(dir));
    
    % 3D-SIM: kz1
    for index = 0:1
        z_top_1 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
        z_bottom_1 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz1;
        h_top = fsurf(x_negative, y_negative, z_top_1, [theta_min, theta_max, phi_min, phi_max]);
        h_bottom = fsurf(x_negative, y_negative, z_bottom_1, [theta_min, theta_max, phi_min, phi_max]);
        
        h_top.FaceColor = [0 0.4470 0.7410];
        h_top.FaceAlpha = 1;
        h_top.EdgeColor = 'none';
        h_top.AmbientStrength = 0.3;
        
        h_bottom.FaceColor = [0 0.4470 0.7410];
        h_bottom.FaceAlpha = 1;
        h_bottom.EdgeColor = 'none';
        h_bottom.AmbientStrength = 0.3;
    end
    
    % 3D-SIM: kz3
    for index = 0:1
        z_top_1 = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz3;
        z_bottom_1 = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max)) + (-1)^index * kz3;
        h_top = fsurf(x_negative, y_negative, z_top_1, [theta_min, theta_max, phi_min, phi_max]);
        h_bottom = fsurf(x_negative, y_negative, z_bottom_1, [theta_min, theta_max, phi_min, phi_max]);
        
        h_top.FaceColor = [0 0.4470 0.7410];
        h_top.FaceAlpha = 1;
        h_top.EdgeColor = 'none';
        h_top.AmbientStrength = 0.3;
        
        h_bottom.FaceColor = [0 0.4470 0.7410];
        h_bottom.FaceAlpha = 1;
        h_bottom.EdgeColor = 'none';
        h_bottom.AmbientStrength = 0.3;
    end

    
    %% order = 2
    z_top = @(theta, phi) (RI / lambda_Em) * (cos(phi) - cos(phi_max));
    z_bottom = @(theta, phi) - (RI / lambda_Em) * (cos(phi) - cos(phi_max));
    % positive direction
    x_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) + 2 * kr * cos(angle(dir));
    y_positive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) + 2 * kr * sin(angle(dir));
    h_top = fsurf(x_positive, y_positive, z_top,[theta_min, theta_max, phi_min, phi_max]);
    h_bottom = fsurf(x_positive, y_positive, z_bottom,[theta_min, theta_max, phi_min, phi_max]);
    
    h_top.FaceColor = [0 0.4470 0.7410];
    h_top.FaceAlpha = 1;
    h_top.EdgeColor = 'none';
    h_top.AmbientStrength = 0.3;
    
    h_bottom.FaceColor = [0 0.4470 0.7410];
    h_bottom.FaceAlpha = 1;
    h_bottom.EdgeColor = 'none';
    h_bottom.AmbientStrength = 0.3;
    
    % negative direction
    x_negetive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * cos(theta) - 2 * kr * cos(angle(dir));
    y_negetive = @(theta, phi) (NA - RI * sin(phi)) / lambda_Em * sin(theta) - 2 * kr * sin(angle(dir));
    h_top = fsurf(x_negetive, y_negetive, z_top,[theta_min, theta_max, phi_min, phi_max]);
    h_bottom = fsurf(x_negetive, y_negetive, z_bottom,[theta_min, theta_max, phi_min, phi_max]);
    
    h_top.FaceColor = [0 0.4470 0.7410];
    h_top.FaceAlpha = 1;
    h_top.EdgeColor = 'none';
    h_top.AmbientStrength = 0.3;
    
    h_bottom.FaceColor = [0 0.4470 0.7410];
    h_bottom.FaceAlpha = 1;
    h_bottom.EdgeColor = 'none';
    h_bottom.AmbientStrength = 0.3;    
end
hold off;

axis equal;
axis off;
axis([-k_max, k_max, -k_max, k_max, -k_max, k_max]);
camlight('right');



