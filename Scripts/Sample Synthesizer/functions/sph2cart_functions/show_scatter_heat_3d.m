% 3d heat scatter
function show_scatter_heat_3d(cam_rft,ct_coord,radar_heat,max_intensity,show_car,th1,th2,th3)
    threshold = max_intensity;
    threshold1 = threshold/th1;
    threshold2 = threshold/th2;
    threshold3 = threshold/th3;
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
    idx_selected_pt = find(radar_heat>=threshold1);
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
    
    % phase 2
    idx_selected_pt = find(radar_heat>=threshold2);
    selected_pt = zeros(size(idx_selected_pt,1),3);
    markerColors_select = zeros(size(idx_selected_pt,1),3);
    for i = 1:size(idx_selected_pt,1)
        selected_pt(i,:) = ct_coord(idx_selected_pt(i),:);
        markerColors_select(i,:) = markerColors(idx_selected_pt(i,:),:);
    end
    x_sl = selected_pt(:,1); y_sl = selected_pt(:,2); z_sl = selected_pt(:,3);
    %ct_coord_ds = ct_coord(1:2000:end,:);
    hold on;
    %scatter3(x_ct,y_ct,z_ct,3,markerColors); hold on;
    scatter3(x_sl,y_sl,z_sl,2,markerColors_select,'filled'); %colorbar();
    
    % phase 3
    idx_selected_pt = find(radar_heat>=threshold3);
    selected_pt = zeros(size(idx_selected_pt,1),3);
    markerColors_select = zeros(size(idx_selected_pt,1),3);
    for i = 1:size(idx_selected_pt,1)
        selected_pt(i,:) = ct_coord(idx_selected_pt(i),:);
        markerColors_select(i,:) = markerColors(idx_selected_pt(i,:),:);
    end
    x_sl = selected_pt(:,1); y_sl = selected_pt(:,2); z_sl = selected_pt(:,3);
    %ct_coord_ds = ct_coord(1:2000:end,:);
    hold on;
    %scatter3(x_ct,y_ct,z_ct,3,markerColors); hold on;
    scatter3(x_sl,y_sl,z_sl,2,markerColors_select,'filled'); %colorbar();
    
    
    colormap(jet); c = colorbar; %caxis([0 1e11]);
    max_tick = strcat(num2str(round(max_intensity/1e12)),'e12');
    c.Ticks = [0 1]; c.TickLabels = {'0',max_tick};
    if show_car == 1
        hold on; 
        %scatter3(ct_coord_ds(:,1),ct_coord_ds(:,2),ct_coord_ds(:,3),0.1,'w'); hold on;
        scatter3(cam_rft(:,1),cam_rft(:,2),cam_rft(:,3),0.5,'filled','k');
    end
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-5 5]); ylim([0 10]); zlim([-1.25,3]); set(gca,'FontSize',font_size);
    view(15,15)
end