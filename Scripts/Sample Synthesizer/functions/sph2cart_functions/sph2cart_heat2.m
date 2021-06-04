% match the heatmap into cartesian coordinates
function heatmap_ct = sph2cart_heat2(scene_lim,N_x,N_y,N_z,grid_pts,grid_heat)
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
    clss = knnsearch(grid_centers,grid_pts,'K',1); % classification
    local_stat = @(x)max(x); % defintion of local statistic
    %class_stat = accumarray(clss,radar_heat,[numr*numc*256 1],local_stat); % data_grouping
    class_stat = accumarray(clss,grid_heat,[N_x*N_y*N_z 1],local_stat); % data_grouping
    heatmap_ct  = reshape(class_stat , size(X)); % 3D reshaping
end