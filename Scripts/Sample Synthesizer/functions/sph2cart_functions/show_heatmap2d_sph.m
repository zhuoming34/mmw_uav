% plot 2d radar heatmaps in spherical coordinates
function show_heatmap2d_sph(heatmap,cam_rft)
    maxheat = max(max(max(heatmap)));
    figure(); 
    font_size = 8;
    
    % Visulize the camera reflectors
    subplot(221); scatter3(cam_rft(:,1),cam_rft(:,2),cam_rft(:,3),0.5,'filled','k');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-4 4]); ylim([0 10]); set(gca,'FontSize',font_size);      
    view(15,30)
    
    xt = linspace(1,64,5);  % az
    %xticks(xt); xticklabels({'30','60','90','120','150'})
    yt = linspace(1,256,11); % rg
    %yticks(yt); yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    zt = linspace(1,64,7);  % el
    %zticks(zt); zticklabels({'75','80','85','90','95','100','105'});

    % Visulize the radar heatmap side view
    radar_heatmap_side = squeeze(max(heatmap,[],2));
    subplot(222); imagesc(radar_heatmap_side.');
    title("Side View"), %set(gca,'XDir','reverse');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('Range'); ylabel('Elevation'); set(gca,'FontSize',font_size);
    xticks(yt); xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    yticks(zt); yticklabels({'75','80','85','90','95','100','105'});
    
    % Visulize the radar heatmap front view
    radar_heatmap_front = squeeze(max(heatmap,[],1));
    subplot(223); imagesc(radar_heatmap_front.');
    title("Front View"), set(gca,'XDir','reverse');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('Azimuth'); ylabel('Elevation'); set(gca,'FontSize',font_size);
    xticks(xt); xticklabels({'30','60','90','120','150'})
    yticks(zt); yticklabels({'75','80','85','90','95','100','105'});
    
    % Visulize the radar heatmap top view
    radar_heatmap_top = squeeze(max(heatmap,[],3));
    subplot(224); imagesc(radar_heatmap_top);
    title("Top View"), set(gca,'XDir','reverse'); set(gca,'YDir','normal');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('Azimuth'); ylabel('Range'); set(gca,'FontSize',font_size);
    xticks(xt); xticklabels({'30','60','90','120','150'})
    yticks(yt); yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
end