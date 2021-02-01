function [] = main
% 2021/01/09 
% new depth settings as zed mini image output HD720
%function [radar_heatmap, visible_cart_v] = main
    % Copyright (c) 2018-2020 Junfeng Guan, Sohrab Madani, Suraj Jog, Saurabh Gupta, 
    % Haitham Hassanieh, University of Illinois at Urbana-Champaign
    % 
    % Permission is hereby granted, free of charge, to any person obtaining a copy
    % of this software and associated documentation files (the "Software"), to deal
    % in the Software without restriction, including without limitation the rights
    % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    % copies of the Software, and to permit persons to whom the Software is
    % furnished to do so, subject to the following conditions:
    % 
    % The above copyright notice and this permission notice shall be included in
    % all copies or substantial port ions of the Software.
    % 
    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    % THE SOFTWARE.
    
    close all; clear; clc;
    
    format shortg; % e.g. 2021 1 10 21 16 28, = Jan 10th, 9:16:28pm 
    clk0 = clock;
    disp("Start");
    disp(clk0);
    
    addpath('functions');
    variable_library;
    
    % depth map resolution, ZED mini HD720 image output -> 1280*720(16:9)
    % video output 2560*720, 0.1-15m, max90*60*100(HxVxD)
    % !but! for HD720 images, HFOV=85deg, VFOV=54deg 
    % L=distance to object(center) in y, image size = (2L)*(2L/sqrt3) 
    numr = 720; numc = 1280;
    os = 5; %1000%2000; % index offset
    %figaddr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/figs/';
    %mataddr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/mats/';
    figaddr = 'F:\3_Education\UMASS\Courses\droneSLAM\mmWave\hawkeye synthesizer\figs\';
    mataddr = 'F:\3_Education\UMASS\Courses\droneSLAM\mmWave\hawkeye synthesizer\mats\';
    
    less_fig = 0;
    show_plot = 0;
    if show_plot == 1
        sp = 'on';
    else
        sp = 'off';
    end
%%
    for CAD_idx = 1:N_CAD_car
        
        % load the surface model
        load(sprintf('../CAD/CAD_model_%d.mat',CAD_idx));
        %load(sprintf('../hawkeye synthesizer/GenerateCADs/CAD_box_%d.mat',CAD_idx));
        %cart_v = myBox;
        %load(sprintf('../hawkeye synthesizer/GenerateCADs/CAD_2boxes_%d.mat',CAD_idx));
        %cart_v = TwoBoxes;
        
        % CAD models are loaded as point clouds of size N_pt by 3, where N_pt
        % is the number of points and 3 values are the cartesian coordinates
        % unit is mm
        
        % Visulize the original point cloud
%         figure; 
%         cart_v_plot = cart_v;
% %         cart_v_plot = datasample(cart_v, 1000); % downsampling when plotting
%         scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
%         xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); axis equal;
%         set(gca,'FontSize',30) % Creates an axes and sets its FontSize to 18

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

        for ks = 1:N_placement_car
            car_scene_v = car1_v_origin;

            %% Rotate     
            car_scene_v.rotate = rotate_ang(randi(length(rotate_ang))); % randomly select a rotation angle and store it in the pc structure
            car_scene_v.rotate = mod(car_scene_v.rotate*(randi(1)*2-1),180);

            
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
            
            translate_x = translate_x_rng(randi(length(translate_x_rng))); % randomly select a translation distance along x axis
            translate_y = translate_y_rng(randi(length(translate_y_rng))); % randomly select a translation distance along y axis
            translate_z = -1250; % translate the point cloud -1250mm to compensate for the height of our radar 

            % translate
            car_scene_v.translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
            car_scene_v.cart_v = car_scene_v.cart_v + car_scene_v.translate; % translate the point cloud
            car_scene_v.bbox = car_scene_v.bbox + car_scene_v.translate; % translate the bounding box
            car_scene_v.lim = [min(car_scene_v.cart_v);max(car_scene_v.cart_v)]; % update the limits in all three dimensions
            
            % convert unit from mm to m
            car_scene_v.cart_v = car_scene_v.cart_v/1000; 
            car_scene_v.bbox = car_scene_v.bbox/1000; 

%             % Visulize the rotated and translated point cloud
            
            cart_v_plot = car_scene_v.cart_v; % downsampling when plotting
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_translatePointCloud','.mat'), 'cart_v_plot');
            if less_fig == 0
                f1 = figure('visible', sp);
                scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
                scatter3(car_scene_v.bbox(:,1), car_scene_v.bbox(:,2),car_scene_v.bbox(:,3),'r'); hold on;
                xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
                xlim([-4 4]); ylim([0 12]);
                set(gca,'FontSize',16); % Creates an axes and sets its FontSize to 18
                saveas(f1, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_1PointCloud"), 'png');
                %savefig(f1,strcat(svaddr,num2str(CAD_idx),"_pm_",num2str(ks+os),"_1PointCloud")); % be careful with the size
            end
            
            
             %% new dep image, 01/09
            [visible_cart_v_dep] = remove_occlusion_radAcam(car_scene_v,"cam",1); % remove occluded body of the car for dep image
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_CameraReflector','.mat'), 'visible_cart_v_dep');
            
            if less_fig == 0
                f2 = figure('visible', sp); 
                cart_v_plot = visible_cart_v_dep; % downsampling when plotting
                scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
                xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
                xlim([-4 4]); ylim([0 12]);
                set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 16
                saveas(f2, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_2DepReflector"), 'png');
            end
            
            ptCloud = visible_cart_v_dep;
            
            pty0 = ptCloud(:,2);
            f = (max(pty0)+min(pty0))/2; % focal lenght in meter
            H_deg = 85; V_deg = 54; % degree, with HD720 
            h = f*tan(H_deg/2/180*pi); v = f*tan(V_deg/2/180*pi); % in meter
            lmin = 0.1; lmax = 15.0; % min/max range of zed mini, in meter
            
            % remove pts that out of bound
            if (min(ptCloud(:,1))<-h || max(ptCloud(:,1))>h ...
                || min(ptCloud(:,3))<-v || max(ptCloud(:,3))>v ...
                || min(ptCloud(:,2))<lmin || max(ptCloud(:,2))>lmax)
                disp("remove pts that out of bound");
                %disp(size(ptCloud2))
                pt_temp = zeros(size(ptCloud,1),1);
                for i = 1:size(ptCloud,1)
                    if ((ptCloud(i,1)<-h)||(ptCloud(i,1)>h) ...
                        || (ptCloud(i,3)<-v)||(ptCloud(i,3)>v) ...
                        || (ptCloud(i,2)<lmin)||(ptCloud(i,2)>lmax))
                        continue;
                    end
                    pt_temp(i) = 1;
                end
                visible_ptCloud_idx = find(pt_temp);
                visible_ptCloud = zeros(size(visible_ptCloud_idx,1),3);
                for i = 1:size(visible_ptCloud_idx,1)
                    visible_ptCloud(i) = ptCloud(visible_ptCloud_idx(i));
                end
            else
                visible_ptCloud = ptCloud;
            end
            
            % corner points 
            ptCorner = [-h,f,-v; -h,f,v; h,f,-v; h,f,v];
            ptCloud2 = [visible_ptCloud; ptCorner];
            
            ptx = ptCloud2(:,1);
            pty = ptCloud2(:,2);
            ptz = ptCloud2(:,3);

            dep0 = pc2i_v2(ptx,pty,ptz,numr,numc,lmin,lmax);
            dep = abs(dep0 - 1)*255; % color flipped
            %disp("got depth info");
            %f3 = figure('visible', sp);
            depcmap = gray; % jet
            %imshow(Dep,depcmap); grid off;
            %colorbar; %colormap gray; %colormap(f3,flipud(gray));       
            %outputBaseFileName = strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_3depthMap");
            outputBaseFileName = strcat(figaddr,num2str(ks+os),'.png');
            imwrite(dep, depcmap, outputBaseFileName);
            %saveas(f3, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_3depthMap"), 'png');
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_depthmap','.mat'), 'dep');
            %disp("saved depth image");
            
            %% Modle radar point reflectors in the scene
            %[visible_cart_v] = remove_occlusion(car_scene_v); % remove occluded body of the car
            [visible_cart_v_rad] = remove_occlusion_radAcam(car_scene_v,"rad",0); 
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_RadarVisible','.mat'), 'visible_cart_v_rad');
            % reflector for radar signal
            try
                reflector_cart_v = model_point_reflector(visible_cart_v_rad,car_scene_v.bbox); % model point reflectors that reflect back to the radar receiver
            catch
                continue;
            end
            if isempty(reflector_cart_v)
                continue;
            end
            % Visulize the radar point reflectors
            cart_v_plot = reflector_cart_v; % downsampling when plotting
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_RadarReflector','.mat'), 'reflector_cart_v');
            if less_fig == 0
                f4 = figure('visible', sp);
                scatter3(cart_v_plot(:,1),cart_v_plot(:,2),cart_v_plot(:,3),10,'filled','k'); hold on;
                xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
                xlim([-4 4]); ylim([0 12]);
                set(gca,'FontSize',16); % Creates an axes and sets its FontSize to 18
                %svaddr = 'F:\3_Education\UMASS\Courses\droneSLAM\mmWave\hawkeye synthesizer\figs\md_';
                saveas(f4, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_4RadarPointReflectors"), 'png');
                %savefig(f2,strcat(svaddr,num2str(CAD_idx),"_pm_",num2str(ks),"_2RadarPointReflectors"));
            end
            %% Simualte received radar signal in the receiver antenna array            
            signal_array = simulate_radar_signal(reflector_cart_v);
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_signal_array','.mat'), 'signal_array');
            %% Radar signal processing, generating 3D radar heatmaps
            radar_heatmap = radar_dsp(signal_array);
            save(strcat(mataddr,'md_',num2str(CAD_idx),'_pm_',num2str(ks+os),'_radar_heatmap','.mat'), 'radar_heatmap');
            %{
            %}
            if less_fig == 0
               % Visulize the radar heatmap top view
                radar_heatmap_top = squeeze(max(radar_heatmap,[],3));
                f5 = figure('visible', sp);
                imagesc(radar_heatmap_top);    
                set(gca,'XDir','reverse');
                set(gca,'YDir','normal');
                colormap jet; %caxis([0 1e11]);
                colorbar;
                %xlabel('Range'); ylabel('Azimuth');
                xlabel('Azimuth'); ylabel('Range'); % 20201129
                set(gca,'FontSize',16); % Creates an axes and sets its FontSize to 18
                saveas(f5, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_5RadarHeatmapTopView"), 'png');
            
                % Visulize the radar heatmap front view
                radar_heatmap_front = squeeze(max(radar_heatmap,[],1));
                f6 = figure('visible', sp);
                imagesc(radar_heatmap_front.');    
                set(gca,'XDir','reverse');
                colormap jet; %caxis([0 1e11]);
                colorbar;
                xlabel('Azimuth'); ylabel('Elevation');
                set(gca,'FontSize',16); % Creates an axes and sets its FontSize to 18
                saveas(f6, strcat(figaddr,'md_',num2str(CAD_idx),"_pm_",num2str(ks+os),"_6RadarHeatmapFrontView"), 'png'); 
            end
            disp(strcat("Model ", num2str(CAD_idx),", placement ", num2str(ks+os), " finished"));
            %format shortg
            clk = clock;
            disp(clk);
        end
    end 
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

