function [] = main
    % 2021/05/25
    % Generate 4-view data including 2D depth images and 3D intensity maps
    % based on Hawkeye radar synthesizer
    % integrated all fucntions including depth images geneartion, heatmap conversion. 
    close all; clear; clc;
    
    format shortg; % e.g. 2021 1 10 21 16 28, = Jan 10th, 9:16:28pm 
    clk0 = clock;disp("Start");disp(clk0);
    
    addpath('functions');
    variable_library;
        
    cam_x = [0,5000,0,-5000]; %mm
    cam_y = [0,5000,10000,5000]; %mm
    cam_ang_deg = [0,90,180,270];
    
    os = 0; %1000%2000; % index offset
    
    rootaddr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/';

    show_plot = 0; sp = 'off';
    if show_plot == 1 
        sp = 'on';
    end
    
    %%% heatmap spherical-to-Cartesian conversion preparation 
    N_phi = 64; N_rho = 256; N_theta = 64;
    scene_lim = [-4, 4; 3, 10; -1.25, 2.75];
    % convert points into cartesian coordinates
    [x_ct,y_ct,z_ct] = sph2cart_pts(N_phi,N_rho,N_theta);
    ct_coord = [x_ct,y_ct,z_ct];

%%
    for CAD_idx = 7
        
        % Routes for saving intermediate/final products
        rftaddr = strcat(rootaddr,'model_',num2str(CAD_idx),'/reflector/');
        sigaddr = strcat(rootaddr,'model_',num2str(CAD_idx),'/signal/');
        figaddr = strcat(rootaddr,'model_',num2str(CAD_idx),'/fig/');
        heataddr = strcat(rootaddr,'model_',num2str(CAD_idx),'/heat2ss/'); % 2ss = 2-snapshot
        cartaddr = strcat(rootaddr,'model_',num2str(CAD_idx),'/cart2ss/'); % converted Cartesian
    
        % load the surface model
        load(sprintf('../CAD/CAD_model_%d.mat',CAD_idx));
      
        % CAD models are loaded as point clouds of size N_pt by 3, where N_pt
        % is the number of points and 3 values are the cartesian coordinates unit is mm
%{        
        % Visulize the original point cloud
        figure; 
        cart_v_plot = cart_v;
        scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),0.5,'filled','k'); hold on;
        xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); axis equal;
        set(gca,'FontSize',30) % Creates an axes and sets its FontSize to 18
%}
        % store point cloud in pc (point cloud) structure
        car_v = car_v_struct;
        car_v.CAD_idx = CAD_idx;
        car_v.N_pt = length(cart_v);
        car_v.cart_v = cart_v;
        car_v.lim = [min(cart_v);max(cart_v)]; % find the limits in all three dimensions 
        [bbox_x, bbox_y, bbox_z] = meshgrid(car_v.lim(:,1),car_v.lim(:,2),car_v.lim(:,3)); % 8 vertices of the bounding box of the point cloud
        car_v.bbox = [bbox_x(:), bbox_y(:), bbox_z(:)]; 
        clear cart_v bbox N_pt car_idx;
        car1_v_origin = car_v;

        for kso = 1:N_placement_car
            ks = kso;
            for cs = 1:4
                car_scene_v = car1_v_origin;
                %% Rotate
                if cs == 1
                    new_rotate = rotate_ang(randi(length(rotate_ang))); % randomly select a rotation angle and store it in the pc structure
                    car_scene_v.rotate = new_rotate;
                     %car_scene_v.rotate = mod(car_scene_v.rotate*(randi(1)*2-1),180);
                else
                    car_scene_v.rotate = new_rotate + cam_ang_deg(cs);
                end
                
                % inline function for 2D rotation
                rotate2d =  @(x, M) (x(:, 1:2) * M);
                rotate_angle_rad = car_scene_v.rotate/180*pi;
                rotation_matrix = [cos(rotate_angle_rad), -sin(rotate_angle_rad); sin(rotate_angle_rad), cos(rotate_angle_rad)]; % create rotation matrix

                car_scene_v.cart_v(:,1:2) = rotate2d(car_scene_v.cart_v, rotation_matrix); % rotate the point cloud 
                car_scene_v.bbox(:,1:2) = rotate2d(car_scene_v.bbox, rotation_matrix); % rotate the bounding box
                car_scene_v.lim = [min(car_scene_v.cart_v);max(car_scene_v.cart_v)]; % update the limits in all three dimensions
                %% Translation
                translate_x_rng = (translate_lim(1,1) - car_scene_v.lim(1,1)):translate_x_res:(translate_lim(1,2) - car_scene_v.lim(2,1)); % range of translation along x axis
                translate_y_rng = (translate_lim(2,1) - car_scene_v.lim(1,2)):translate_y_res:(translate_lim(2,2) - car_scene_v.lim(2,2)); % range of translation along y axis
                %disp('trans lim')disp(translate_lim)disp('car scene lim') disp(car_scene_v.lim)disp('trans x rng')disp(translate_x_rng)disp('trans y rng')disp(translate_y_rng)
                switch cs 
                    case 2               
                        translate_x = new_y - cam_y(cs);
                        translate_y = cam_x(cs) - new_x;
                    case 3               
                        translate_x = -new_x;
                        translate_y = cam_y(cs) - new_y;
                    case 4               
                        translate_x = cam_y(cs) - new_y;
                        translate_y = -cam_x(cs) + new_x;
                    otherwise
                        % case 1
                        new_x = translate_x_rng(randi(length(translate_x_rng))); % randomly select a translation distance along x axis
                        new_y = translate_y_rng(randi(length(translate_y_rng))); % randomly select a translation distance along y axis
                        translate_x = new_x;
                        translate_y = new_y;
                end
                translate_z = -1250; % translate the point cloud -1250mm to compensate for the height of our radar 

                % translate
                car_scene_v.translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
                car_scene_v.cart_v = car_scene_v.cart_v + car_scene_v.translate; % translate the point cloud
                car_scene_v.bbox = car_scene_v.bbox + car_scene_v.translate; % translate the bounding box
                car_scene_v.lim = [min(car_scene_v.cart_v);max(car_scene_v.cart_v)]; % update the limits in all three dimensions

                % convert unit from mm to m
                car_scene_v.cart_v = car_scene_v.cart_v/1000; 
                car_scene_v.bbox = car_scene_v.bbox/1000; 

                %%% ---------- Depth Image ----------
                % select camera reflectors
                [visible_cart_v_dep] = remove_occlusion_v2(car_scene_v,"cam",0); % remove occluded body of the car for dep image
                save(strcat(rftaddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),"_cam_",num2str(cs),'_CameraReflector','.mat'), 'visible_cart_v_dep');
                [Dep1,Dep2,reszImg] = pc2dep(visible_cart_v_dep); % generate depth images
                % save depth images
                cmap1 = jet; cmap2 = gray;
                outputImgName1 = strcat(figaddr,'original/cam',num2str(cam),'/',num2str(idx+os),'.png');      
                outputImgName2 = strcat(figaddr,'extend/cam',num2str(cam),'/',num2str(idx+os),'.png');            
                outputImgName3 = strcat(figaddr,'color128/cam',num2str(cam),'/',num2str(idx+os),'.png');
                outputImgName4 = strcat(figaddr,'gray128/cam',num2str(cam),'/',num2str(idx+os),'.png');
                imwrite(Dep1, cmap1, outputImgName1); 
                imwrite(Dep2, cmap1, outputImgName2); 
                imwrite(reszImg, cmap1, outputImgName3);
                imwrite(reszImg, cmap2, outputImgName4); 
                
                %%% ---------- Radar Heatmap ----------
                % Modle radar point reflectors in the scene
                %[visible_cart_v] = remove_occlusion(car_scene_v); % remove occluded body of the car
                [visible_cart_v_rad] = remove_occlusion_v2(car_scene_v,"rad",0); 
                try
                    reflector_cart_v = model_point_reflector(visible_cart_v_rad,car_scene_v.bbox); % model point reflectors that reflect back to the radar receiver
                catch
                    continue; % disp('error'); 
                end
                if isempty(reflector_cart_v)
                    continue;
                end

                % adding environmental noise
                [evnoise] = add_evn_noise();
                reflector_cart_v_noisy = [reflector_cart_v;evnoise];
                save(strcat(rftaddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_cam_',num2str(cs),'_RadarReflectorNoisy','.mat'), 'reflector_cart_v_noisy');

                % Simualte received radar signal in the receiver antenna array     
                % signal_array = simulate_radar_signal(reflector_cart_v); % with environment noise
                signal_array_noisy = simulate_radar_signal(reflector_cart_v_noisy); % with environment noise
                %save(strcat(sigaddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_cam_',num2str(cs),'_signal_array_Noisy','.mat'), 'signal_array_noisy');
                
                % Radar signal processing, generating 3D radar heatmaps
                % radar_heatmap_noisy = radar_dsp(signal_array_noisy); % full-scale
                signal_array_noisy2 = signal_array_noisy(:,:,32:33); % 2 snapshots in the middle vertically
                radar_heatmap_noisy = radar_dsp2ss(signal_array_noisy2);
                %save(strcat(heataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_cam_',num2str(cs),'_radar_heatmap2_noisy','.mat'), 'radar_heatmap_noisy');

                % Convert spherical coordinate to Cartesian coordinate
                % with a specified boundary
                heatmap_ct = heatSph2Cart(radar_heatmap_noisy, ct_coord, scene_lim, N_phi, N_rho, N_theta);
                save(strcat(cartaddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_cam_',num2str(cs),'_cart_heatmap2_noisy','.mat'), 'heatmap_ct');
                
                disp(strcat("Model ", num2str(CAD_idx),", placement ", num2str(ks+os),", camera ", num2str(cs), " finished"));           
                clk = clock; disp(clk); %format shortg
            end
        end
    end 
    
    show_time(clk0)
end

function show_time(clk0)
    disp("started from");
    disp(clk0);
    disp("All finished, total used");
    clk1 = clock-clk0;
    if clk1(6)<0 
        clk1(6)=clk1(6)+60;
        clk1(5) = clk1(5) - 1;
    end
    if clk1(5)<0
        clk1(5)=clk1(5)+60;
        clk1(4) = clk1(4) - 1;  
    end
    if clk1(4)<0
        clk1(4)=clk1(4)+24;
        % back-counting by hrs
    end
    disp(clk1);
end
