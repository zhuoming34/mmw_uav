% 01/17/2021
% ver2.3, for converting multiple .mat files
% use pts with intensities that are larger than 1/10 of the largest
% intensity only, and limit the size of meshgrid to
% x=[-4,4],y=[3,10],z=[-1.25,2.75]
% convert radar heatmap in hawkeye synthesizer
% from spherical coordinates (r,theta,phi)to cartesian coordinates (x,y,z)
close all; clear; clc;

N_phi = 64; N_rho = 256; N_theta = 64; 
% convert points into cartesian coordinates
[x_ct,y_ct,z_ct] = sph2cart_pts(N_phi,N_rho,N_theta);
ct_coord = [x_ct,y_ct,z_ct];

addr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/20210509/heat2ss/';
%addr = 'F:/3_Education/UMASS/Courses/droneSLAM/mmWave/hawkeye synthesizer/dataset/20210103/';
saveaddr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/20210509/cart2ss2/';
CAD_idx = 0;
offset = 0;
camos = 0;
for idx_mat = 1:600
    for cam = 1:4
        filename = strcat('md_',num2str(CAD_idx),'_pm_',num2str(idx_mat+offset),'_cam_',num2str(cam+camos),'_radar_heatmap2_noisy');
        mat = load(strcat(addr,filename,'.mat'));
        heatmap = mat.radar_heatmap_noisy; % rho*azimuth*elevation
        % match intensity values to corresponding points
        radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta);

        % gets points whose intensities are larger than the threshold
        threshold = max(max(max(heatmap)))/10;
        idx_heat_fl = find(radar_heat >= threshold);
        pt_heat_fl = zeros(size(idx_heat_fl,1),3);
        heat_fl = zeros(size(idx_heat_fl,1),1);
        for i = 1:size(idx_heat_fl,1)
            pt_heat_fl(i,:) = ct_coord(idx_heat_fl(i),:);
            heat_fl(i) = radar_heat(idx_heat_fl(i));
        end

        % no filtering
        %pt_heat_fl = ct_coord;
        %heat_fl = radar_heat;
        %x_hf = pt_heat_fl(:,1); y_hf = pt_heat_fl(:,2); z_hf = pt_heat_fl(:,3);
        
        % convert heatmap from spherical to cartesian coordinate system
        N_x= 64; N_y = 256; N_z = 64;
        %scene_lim = [-4, 4; 3, 10; -1.25, 2.75];
        %scene_lim = [-2.5, 2.5; 0, 5; -0.5, 2.5]; % boxes
        scene_lim = [-1, 1; 1.5, 5; -0.5, 1]; % boxes 2
        heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,pt_heat_fl,heat_fl);
        save(strcat(saveaddr,'md_',num2str(CAD_idx),'_pm_',num2str(idx_mat+offset),'_cam_',num2str(cam+camos),'_cart_heatmap2_noisy','.mat'), 'heatmap_ct');
        disp(strcat("Model ", num2str(CAD_idx),", placement ", num2str(idx_mat+offset),', camera ',num2str(cam+camos), ": done"));
    end
end
disp('finished');
%% plots
show_plots = 0; idx_plots = [4];
if show_plots == 1
    for idx_plot = idx_plots
        if idx_plot == 1 || idx_plot == 3
            idx = 8;
            filename2 = strcat('md_1_pm_',num2str(idx),'_CameraReflector');
            mat2 = load(strcat(addr,filename2,'.mat'));
            cam_rft = mat2.visible_cart_v_dep; % camera reflectors
        end
        % plot 2d radar heatmaps in spherical coordinates
        if idx_plot == 1 
            show_heatmap2d(heatmap, cam_rft);
        end
        % plot point cloud in cartesian coordinate system
        if idx_plot == 2
            show_ct_pt(ct_coord);
        end
        % 3d heat scatter
        if idx_plot == 3
            max_intensity = max(max(max(heatmap)));
            show_car = 0;
            show_scatter_heat_3d(cam_rft,ct_coord,radar_heat,max_intensity,show_car);
        end
        % plot 3d heatmap slice
        if idx_plot == 4
            show_slice_heat_3d(heatmap_ct)
            view(0,0) % front(0,0), top(0,90), left(-90,0),right(90,0)
        end
    end
end
%% functions 

% convert points into cartesian coordinates
function [x_ct,y_ct,z_ct] = sph2cart_pts(N_phi,N_rho,N_theta)
    c = 3e8; % speed of light 
    %BW = 3987.61e6; % Bandwidth = 3987.61MHz
    Fs = 2047e3;
    sweep_slope = 29.982e+12; % Hz/s
    %rho_res = c/2/BW; % in m, range resolution, 3.8cm for 3987.61MHz
    rho_min = 0; rho_max = Fs*c/2/sweep_slope; % range, m
    phi_min = 30; phi_max = 150; % azimuth, deg
    theta_min = 75; theta_max = 105; % elevation, deg
    % get corresponding positions of each cell
    ticks_phi = linspace(phi_min,phi_max,N_phi); % azimuth axis of the output radar heatmap in degree
    ticks_theta = linspace(theta_min,theta_max,N_theta);
    ticks_rho = linspace(rho_min,rho_max,N_rho);
    ticks_phi = ticks_phi/180*pi;
    ticks_theta = ticks_theta/180*pi;
    % to use sph2cart, theta is the elevation angle from x-y plane
    ticks_theta_top = pi/2 - ticks_theta(1:32);
    ticks_theta_bottom = -(ticks_theta(33:64) - pi/2);
    ticks_theta = [ticks_theta_top,ticks_theta_bottom];

    N_cell = N_phi*N_theta*N_rho;
    sp_coord = zeros(N_cell,3);
    idx_sp_coord = 1;
    for idx_rho = 1:N_rho
        for idx_phi = 1:N_phi
            for idx_theta = 1:N_theta
                sp_coord(idx_sp_coord,:) = [ticks_rho(idx_rho),ticks_phi(idx_phi),ticks_theta(idx_theta)];
                idx_sp_coord = idx_sp_coord + 1;
            end
        end
    end
    [x_ct,y_ct,z_ct] = sph2cart(sp_coord(:,2),sp_coord(:,3),sp_coord(:,1)); % [x,y,z] = sph2cart(azimuth,elevation,r)
end

% match intensity values to corresponding points
function radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta)
    N_cell = N_phi*N_theta*N_rho;
    radar_heat = zeros(N_cell,1);
    idx_sp_coord = 1;
    for idx_rho = 1:N_rho
        for idx_phi = 1:N_phi
            for idx_theta = 1:N_theta
                radar_heat(idx_sp_coord) = heatmap(idx_rho,idx_phi,idx_theta); 
                idx_sp_coord = idx_sp_coord + 1;
            end
        end
    end
end

% match the heatmap into cartesian coordinates
function heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,pts,radar_heat)
    x_min = scene_lim(1,1); x_max = scene_lim(1,2);
    y_min = scene_lim(2,1); y_max = scene_lim(2,2);
    z_min = scene_lim(3,1); z_max = scene_lim(3,2);
    xx = linspace(x_min,x_max,N_x);
    yy = linspace(y_min,y_max,N_y);
    zz = linspace(z_min,z_max,N_z);
    % create a meshgrid and assign points to corresponding cubes
    [X,Y,Z] = meshgrid(xx,yy,zz);
    grid_centers = [X(:),Y(:),Z(:)];
    %clss = knnsearch(grid_centers,[x_ct,y_ct,z_ct]); % classification
    % pts = [x_ct,y_ct,z_ct];
    clss = knnsearch(grid_centers,pts,'K',1); % classification
    local_stat = @(x)mean(x); % defintion of local statistic
    %class_stat = accumarray(clss,radar_heat,[numr*numc*256 1],local_stat); % data_grouping
    class_stat = accumarray(clss,radar_heat,[N_x*N_y*N_z 1],local_stat); % data_grouping
    heatmap_ct  = reshape(class_stat , size(X)); % 3D reshaping

    for id_y = 1:N_y
        for id_x = 1:N_x
            for id_z = 1:N_z
                if heatmap_ct(id_y,id_x,id_z) == 0
                    if id_x == 1
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/2;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/2;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        end
                    elseif id_x == N_x
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/2;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/2;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        end
                    else
                        if id_z == 1
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1))/3;
                        elseif id_z == N_z
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z-1))/3;
                        else
                            heatmap_ct(id_y,id_x,id_z) = (heatmap_ct(id_y,id_x+1,id_z)+heatmap_ct(id_y,id_x-1,id_z)+heatmap_ct(id_y,id_x,id_z+1)+heatmap_ct(id_y,id_x,id_z-1))/4;
                        end
                    end
                end
            end
        end
    end

end

% plot 2d radar heatmaps in spherical coordinates
function show_heatmap2d(heatmap, cam_rft)
    figure(); 
    font_size = 8;
    
    % Visulize the camera reflectors
    subplot(221); scatter3(cam_rft(:,1),cam_rft(:,2),cam_rft(:,3),10,'filled','k');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-4 4]); ylim([0 10]); set(gca,'FontSize',font_size);      
    view(15,30)
    
    % Visulize the radar heatmap side view
    radar_heatmap_side = squeeze(max(heatmap,[],2));
    subplot(222); imagesc(radar_heatmap_side.');
    %set(gca,'XDir','reverse');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Range'); ylabel('Elevation'); set(gca,'FontSize',font_size);
    
    % Visulize the radar heatmap front view
    radar_heatmap_front = squeeze(max(heatmap,[],1));
    subplot(223); imagesc(radar_heatmap_front.');
    set(gca,'XDir','reverse');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Azimuth'); ylabel('Elevation'); set(gca,'FontSize',font_size);
    
    % Visulize the radar heatmap top view
    radar_heatmap_top = squeeze(max(heatmap,[],3));
    subplot(224); imagesc(radar_heatmap_top);
    set(gca,'XDir','reverse'); set(gca,'YDir','normal');
    colormap jet; %caxis([0 1e11]); %colorbar;
    xlabel('Azimuth'); ylabel('Range'); set(gca,'FontSize',font_size);

end

% plot point cloud in cartesian coordinate system
function show_ct_pt(ct_coord)
    font_size = 8;
    figure(); scatter3(ct_coord(:,1),ct_coord(:,2),ct_coord(:,3),1,'filled','k');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-10 10]); ylim([0 10]); zlim([-5,5]); set(gca,'FontSize',font_size);      
    %view(15,30)
end

% plot 3d heatmap slice
function show_slice_heat_3d(heatmap_ct)
    m = size(heatmap_ct, 2); % az
    n = size(heatmap_ct, 3); % el
    l = size(heatmap_ct, 1); % rg
    xi = linspace(1,l,l); % range
    yi = linspace(1,m,m); % azimuth
    zi = linspace(1,n,n); % elevation
    [XX,YY,ZZ] = meshgrid(yi,xi,zi); % [l * m * n]

    xslice = 1:m;    % location of y-z planes
    yslice = 1:l;    % location of x-z plane
    zslice = 1:n;    % location of x-y planes
    
    figure(); h = slice(XX,YY,ZZ,heatmap_ct,xslice,yslice,zslice);
    xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
    xlim([0 64]); ylim([0 256]); zlim([0,64]);
    xt = linspace(0,64,9); xticks(xt); 
    xticklabels({'-4','-3','-2','-1','0','1','2','3','4'})
    yt = linspace(0,256,8); yticks(yt); 
    yticklabels({'3','4','5','6','7','8','9','10'})
    zt = linspace(0,64,17); zticks(zt); 
    zticklabels({'-1.25','-1','-0.75','-0.5','-0.25','-0','0.25',...
        '0.5','0.75','1','1.25','1.5','1.75','2','2.25','2.5','2.75'});
    %set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');
    set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','flat');
    %set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','0.01');
    % set transparency to correlate to the data values.
    alpha('color'); %alpha(h, 0.01);
    colorbar; colormap jet;
end

% 3d heat scatter
function show_scatter_heat_3d(cam_rft,ct_coord,radar_heat,max_intensity,show_car)
    threshold = max_intensity/10;
    cmap = jet;
    %radar_heat_sort = sort(radar_heat);
    v = rescale(radar_heat, 1, 256); % Nifty trick!
    numValues = length(radar_heat);
    markerColors = zeros(numValues, 3);
    % Now assign marker colors according to the value of the data.
    for k = 1 : numValues
        row = round(v(k));
        markerColors(k, :) = cmap(row, :);
    end
    % filter out low-intensity points
    %max_intensity = max(max(max(heatmap)));
    idx_selected_pt = find(radar_heat>=threshold/50);
    selected_pt = zeros(size(idx_selected_pt,1),3);
    markerColors_select = zeros(size(idx_selected_pt,1),3);
    for i = 1:size(idx_selected_pt,1)
        selected_pt(i,:) = ct_coord(idx_selected_pt(i),:);
        markerColors_select(i,:) = markerColors(idx_selected_pt(i,:),:);
    end
    x_sl = selected_pt(:,1); y_sl = selected_pt(:,2); z_sl = selected_pt(:,3);
    %ct_coord_ds = ct_coord(1:2000:end,:);
    figure(); font_size = 8;
    %scatter3(x_ct,y_ct,z_ct,3,markerColors); hold on;
    scatter3(x_sl,y_sl,z_sl,10,markerColors_select,'filled'); %colorbar();
    colormap(jet); c = colorbar; %caxis([0 1e11]);
    max_tick = strcat(num2str(round(max_intensity/1e12)),'e12');
    c.Ticks = [0 1]; c.TickLabels = {'0',max_tick};
    if show_car == 1
        hold on; 
        %scatter3(ct_coord_ds(:,1),ct_coord_ds(:,2),ct_coord_ds(:,3),0.1,'w'); hold on;
        scatter3(cam_rft(:,1),cam_rft(:,2),cam_rft(:,3),10,'filled','k');
    end
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-10 10]); ylim([0 10]); zlim([-5,5]); set(gca,'FontSize',font_size);
    view(15,15)
end