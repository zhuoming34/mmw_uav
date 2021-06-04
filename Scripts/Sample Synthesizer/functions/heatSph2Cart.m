function heatmap_ct = heatSph2Cart(heatmap, ct_coord, scene_lim, N_x, N_y, N_z, ptGrid, ptGrid_heat)
        
    addpath('sph2cart_functions');

    % match intensity values to corresponding spherical voxel center
    radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta);

    % ver 2.3
    % get points whose intensities are larger than the threshold
    threshold = max(max(max(heatmap)))/50;
    idx_heat_fliter = find(radar_heat >= threshold);
    ct_coord_filtered = zeros(size(idx_heat_fliter,1),3);
    radar_heat_filtered = zeros(size(idx_heat_fliter,1),1);
    for i = 1:size(idx_heat_fliter,1)
        ct_coord_filtered(i,:) = ct_coord(idx_heat_fliter(i),:);
        radar_heat_filtered(i) = radar_heat(idx_heat_fliter(i));
    end
    % no filtering
    %pt_heat_fl = ct_coord;
    %heat_fl = radar_heat;
    %x_hf = pt_heat_fl(:,1); y_hf = pt_heat_fl(:,2); z_hf = pt_heat_fl(:,3);
    % convert heatmap from spherical to cartesian coordinate system
    % heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,pt_heat_fl,heat_fl);


    % ver 3.1
    %ptGrid_clss = knnsearch(ct_coord,ptGrid,'K',1); % classification
    ptGrid_clss = knnsearch(ct_coord_filtered,ptGrid,'K',1); % classification
    for i = 1:numPoint
        %ptGrid_heat(i) = radar_heat(ptGrid_clss(i));
        ptGrid_heat(i) = radar_heat_filtered(ptGrid_clss(i));
    end 
    heatmap_ct = sph2cart_heat2(scene_lim,N_x,N_y,N_z,ptGrid,ptGrid_heat);

    %{
    % logarithmic scale
    heatmap_ct = heatmap_ct + 1.0;
    heatmap_ct = log(heatmap_ct);

    % ver 3.3 filtering
    maxheat = max(max(max(heatmap_ct)));
    filter_idx = find(heatmap_ct < maxheat*0.675);
    heatmap_ct(filter_idx) = 0;  
    %}   
end