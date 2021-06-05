function heatmap_ct = heatSph2Cart(heatmap, scene_lim, N_x, N_y, N_z, ptGrid, ptGrid_heat, ptGrid_clss)
        
    addpath('sph2cart_functions');

    % match intensity values to corresponding spherical voxel center
    radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta);
    for i = 1:numPoint
        ptGrid_heat(i) = radar_heat(ptGrid_clss(i));
    end 

    % filtering: get points whose intensities are larger than the threshold
    threshold = max(max(max(heatmap)))/50;
    idx_heat_fliter = find(ptGrid_heat >= threshold);
    ptGrid_filtered = zeros(size(idx_heat_fliter,1),3);
    ptGrid_heat_filtered = zeros(size(idx_heat_fliter,1),1);
    for i = 1:size(idx_heat_fliter,1)
        ptGrid_filtered(i,:) = ptGrid(idx_heat_fliter(i),:);
        ptGrid_heat_filtered(i) = ptGrid_heat(idx_heat_fliter(i));
    end

    % generate 3d heatmap in Cartesian
    heatmap_ct = sph2cart_heat2(scene_lim,N_x,N_y,N_z,ptGrid_filtered,ptGrid_heat_filtered);

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