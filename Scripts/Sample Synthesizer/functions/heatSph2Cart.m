function heatmap_ct = heatSph2Cart(heatmap, ct_coord, scene_lim, N_x, N_y, N_z)
        
        addpath('sph2cart_functions');

        % match intensity values to corresponding points
        radar_heat = matchHeat(heatmap,N_phi,N_rho,N_theta);

        % to reduce the size of heatmap
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
        %N_x= 64; N_y = 256; N_z = 64;
        %scene_lim = [-4, 4; 3, 10; -1.25, 2.75];
        heatmap_ct = sph2cart_heat(scene_lim,N_x,N_y,N_z,pt_heat_fl,heat_fl);
end