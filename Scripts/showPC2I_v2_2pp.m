% 01/21/2021
% generate 2d depth image from camera reflector point cloud
% with perspective projection
close all; clear; clc;
% --------------------------------------
% https://support.stereolabs.com/hc/en-us/articles/360007395634-What-is-the-camera-focal-length-and-field-of-view-
% zed mini image output resolution of HD720
% width = 1280, height = 720, pixel number = 1280*720
% has a focal length of 700 pixels
% pixel size of 0.004mm
% field of view: fov = 2*arctan(pixelNumber/(2*focalLength))*(180/pi)
% VFOV = 54deg, HFOV = 85deg
% "The focus distance is fixed. All objects at distances from 28cm out to infinity will be sharp"
% depth range = [0.1, 15] meters
numc = 1280; numr = 720;
f_px = 700; % focal lenght in pixel
pxSize = 0.004; % mm
f = f_px*pxSize/1000; % 2.8mm -> 0.0028m
ppH = numc*pxSize/1000; ppV = numr*pxSize/1000; % % projection plane size in m
hfov = 2*atan(numc/2/f_px); vfov = 2*atan(numr/2/f_px); % rad
hfov_deg = hfov/pi*180; vfov_deg = vfov/pi*180;
r1 = 0.1; r2 = 15.0; % meter
h = r2*tan(hfov/2); v = r2*tan(vfov/2); % max range the camera covers

% scene boundary
%scene_lim = [-4,4; 3,10; -1.25, 2.75];
scene_corners = [-4,3,-1.25; -4,3,2.75;  4,3,-1.25; 4,3,2.75];
%                  4,10,-1.25;  4,10,2.75; -4,10,2.75; -4,10,-1.25];
%pts = [ptCloud;scene_corners];
sn_x = scene_corners(:,1); sn_y = scene_corners(:,2); sn_z = scene_corners(:,3);

% scene projection
sn_py = linspace(f,f,size(sn_y,1))'; % focal length
sn_ratio = sn_y./sn_py;
sn_px = -sn_x./sn_ratio;
sn_pz = -sn_z./sn_ratio;

os = 0; % index offset
addr = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/mats/';
outaddr1 = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/figs/original/';
outaddr2 = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/figs/extended/';
outaddr3 = '/home/huang/Documents/HawkEye-Data-Code-master/Synthesizer/figs/resized/';

for idx = 1:1500
    filename = strcat('md_1_pm_',num2str(idx),'_CameraReflector.mat');
    path = strcat(addr,filename);
    load(path);
    ptCloud = visible_cart_v_dep;
    x = ptCloud(:,1); y = ptCloud(:,2); z = ptCloud(:,3);

    % perspective projection
    py = linspace(f,f,size(y,1))'; % focal length
    ratio = y./py;
    px = -x./ratio;
    pz = -z./ratio;

    % remove pts that out of bound
    if (min(px)<-ppH||max(px)>ppH || min(pz)<-ppV||max(pz)>ppV)
        disp("remove out-of-bound points");
        ppts = [px,py,pz];
        pt_temp = zeros(size(ppts,1),1);
        for i = 1:size(ppts,1)
            if (px(i)<-ppH||px(i)>ppH || pz(i)<-ppV||pz(i)>ppV)
                continue;
            end
            pt_temp(i) = 1;
        end
        visible_ptCloud_idx = find(pt_temp);
        visible_ptCloud = zeros(size(visible_ptCloud_idx,1),3);
        visible_ptCloud_pp = zeros(size(visible_ptCloud_idx,1),3); % pts on projection plane
        for i = 1:size(visible_ptCloud_idx,1)
            visible_ptCloud(i) = ptCloud(visible_ptCloud_idx(i));
            visible_ptCloud_pp(i) = ppts(visible_ptCloud_idx(i));
        end
        x = visible_ptCloud(:,1); y = visible_ptCloud(:,2); z = visible_ptCloud(:,3);
        px = visible_ptCloud_pp(:,1); py = visible_ptCloud_pp(:,2); pz = visible_ptCloud_pp(:,3);
    else
        disp("no point is out of bound");
    end

    % ---- depth calculation ----
    d = sqrt(x.^2 + y.^2 + z.^2);
    %sn_d = sqrt(sn_x.^2 + sn_y.^2 + sn_z.^2);

    % grid construction
    xl = -ppH/2; xr = ppH/2; zl = -ppV/2; zr = ppV/2;
    xx = linspace(xr,xl,numc); zz = linspace(zl,zr,numr);
    [X,Z] = meshgrid(xx,zz);
    grid_centers = [X(:),Z(:)];

    % classification
    clss = knnsearch(grid_centers,[px,pz]); 
    % defintion of local statistic
    local_stat = @(x)min(x); 
    %local_stat = @(x)mean(x); 
    % data_grouping
    class_stat = accumarray(clss,d,[numr*numc 1],local_stat);
    % 2D reshaping
    class_stat_M  = reshape(class_stat , size(X)); 
    % Force un-filled cells to the brightest color
    % add limits at front and back
    s1 = size(class_stat_M,1);
    s2 = size(class_stat_M,2);
    class_stat_M(1,s2/2) = r1; % 0.1m
    class_stat_M(s1,s2/2) = r2; % 15.0m
    class_stat_M (class_stat_M == 0) = max(max(class_stat_M));
    % flip image horizontally and vertically
    %I0 = class_stat_M(end:-1:1,end:-1:1);
    %I0 = class_stat_M(end:-1:1,:);
    %I0 = class_stat_M(:,end:-1:1); % flip image horizontally
    % normalize pixel values to [0,1]
    I0 = class_stat_M;
    I = ( I0 - min(min(I0)) ) ./ ( max(max(I0)) - min(min(I0)) );
    
    % extend and crop to the same size of scene
    N_col_l = ceil(abs((sn_px(3)-xl)/(pxSize/1000))); % extend
    N_col_r = ceil(abs((sn_px(1)-xr)/(pxSize/1000))); % extend
    N_row_t = floor(abs((sn_pz(1)-zr)/(pxSize/1000))); % crop
    N_row_b = ceil(abs((sn_pz(2)-zl)/(pxSize/1000))); % extend
    N_col = numc + N_col_l + N_col_r; N_row = numr - N_row_t + N_row_b;
    col_l = ones(N_row,N_col_l); col_r = ones(N_row,N_col_r);
    row_b = ones(N_row_b, numc);
    %I2 = [col_l,[I(N_row_t+1:end,:);row_b],col_r];
    I2 = [col_l,[row_b;I(1:end-N_row_t,:)],col_r];

    % saving images
    Dep1 = abs(I - 1)*255;
    Dep2 = abs(I2 - 1)*255;
    reszImg = imresize(Dep2, [256,256]);

    outputBaseFileName1 = strcat(outaddr1,num2str(idx+os),'.png');
    imwrite(I, outputBaseFileName1); 

    outputBaseFileNameq2 = strcat(outaddr2,num2str(idx+os),'.png');
    imwrite(Dep, outputBaseFileName2); 

    outputBaseFileName3 = strcat(outaddr3,num2str(idx+os),'.png');
    imwrite(reszImg, outputBaseFileName3); 
    
    disp(idx+os)
end
disp('finished')
%%
show_scene = 0;
if show_scene == 1
    %plot the point cloud
%step = 1000;
%x1 = x(1:step:end,1);
%y1 = y(1:step:end,1);
%z1 = z(1:step:end,1);
scene_corners2 = [-4, 3,-1.25; -4, 3,2.75; 4, 3,2.75;  4, 3,-1.25;
                  4,10,-1.25;  4,10,2.75; -4,10,2.75; -4,10,-1.25;
                  -4, 3,-1.25; 4, 3,-1.25; 4, 3,2.75; 4,10,2.75;
                  4,10,-1.25;  -4,10,-1.25; -4,10,2.75; -4, 3,2.75;];
sn_x2 = scene_corners2(:,1); sn_y2 = scene_corners2(:,2); sn_z2 = scene_corners2(:,3);
   
hfov_rd = 120/180*pi; vfov_rd = 30/180*pi;
r3 = 10.0; % meter
h2 = r3*tan(hfov_rd/2); v2 = r3*tan(vfov_rd/2); % max range the camera covers
figure(1);
scatter3(x,y,z,0.2,'filled','k');
hold on; scatter3(0,0,0,50,'filled','m');
hold on; scatter3(sn_x,sn_y,sn_z,10,'r'); % bounding box corners
%hold on; plot3(sn_x,sn_y,sn_z,':r'); % bounding box lines 1
hold on; plot3(sn_x2,sn_y2,sn_z2,':r'); % bounding box lines 2
hold on; plot3([0 h],[0 r2],[0 0],'-.g');hold on; plot3([0 -h],[0 r2],[0 0],'-.g')
hold on; plot3([0 0],[0 r2],[0 v],'-.g');hold on; plot3([0 0],[0 r2],[0 -v],'-.g')
hold on; plot3([0 h2],[0 r3],[0 0],'-.b');hold on; plot3([0 -h2],[0 r3],[0 0],'-.b')
hold on; plot3([0 0],[0 r3],[0 v2],'-.b');hold on; plot3([0 0],[0 r3],[0 -v2],'-.b')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
xlim([-7,7]), ylim([0,12]), zlim([-3,3])
%view(90,0)
end