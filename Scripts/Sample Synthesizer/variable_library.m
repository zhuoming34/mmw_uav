N_CAD_car=1; % number of CAD models of cars, max =38
N_placement_car = 1428; % # of placement we create with every selected car/group of cars

% angle of rotation for every model
rotate_ang_coarse = [0:5:360];%[5:5:175]; % [0,90]=head left, [90,180]=head right
rotate_ang_fine = [-5:5];
rotate_ang = [];
for k_rotate_ang_fine = rotate_ang_fine
    rotate_ang = [rotate_ang,rotate_ang_coarse+k_rotate_ang_fine];
end
rotate_ang = sort(rotate_ang);

translate_lim = [-3000, 3000; 2000, 8000]; % limits of the translation along the x and y axis
%translate_lim = [-600, 600; 2300, 3500]; % 4-way boundary
translate_x_res = 10; %500; % resolution of translation along the x axis unit: mm
translate_y_res = 10; %500; % resolution of translation along the y axis unit: mm

% structure of point cloud and information of the car CAD model
car_v_struct = struct('cart_v',[], ... % cartesian coordinates vector
                'sph_v',[], ... % spherical coordinates vector
                'N_pt',0, ... % # of points in the model
                'bbox',[], ... % bounding box of the car 
                'lim',[], ... % min and max xyz coordinates
                'CAD_idx',0, ...
                'rotate',[], ... % degree rotated
                'translate',[]); % distance translated


