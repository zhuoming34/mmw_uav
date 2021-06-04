% plot 2d radar heatmaps in Cartesian coordinates
function show_heatmap2d_cart(heatmap,cam_rft)   
    
    figure(); 
    font_size = 8;
    % Visulize the camera reflectors
    subplot(221); scatter3(cam_rft(:,1),cam_rft(:,2),cam_rft(:,3),0.5,'filled','k');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal
    xlim([-4, 4]); ylim([0, 10]); zlim([-1.25, 2.75]); 
    set(gca,'FontSize',font_size);      
    view(15,30)
       
    xt = linspace(1,64,9); yt = linspace(1,256,8); zt = linspace(1,64,17); 
    
    % Visulize the radar heatmap side view
    radar_heatmap_side = squeeze(max(heatmap,[],2));
    subplot(222); imagesc(radar_heatmap_side.');
    title("Side View")
    set(gca,'XDir','normal');set(gca,'YDir','normal');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('y(m)'); ylabel('z(m)'); set(gca,'FontSize',font_size);
    xticks(yt); xticklabels({'3','4','5','6','7','8','9','10'})
    yticks(zt); 
    yticklabels({'-1.25','-1','-0.75','-0.5','-0.25','-0','0.25',...
        '0.5','0.75','1','1.25','1.5','1.75','2','2.25','2.5','2.75'});
    
    
    % Visulize the radar heatmap front view
    radar_heatmap_front = squeeze(max(heatmap,[],1));
    subplot(223); imagesc(radar_heatmap_front.');
    title("Front View")
    %set(gca,'XDir','reverse');
    set(gca,'XDir','normal');set(gca,'YDir','normal');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('x(m)'); ylabel('z(m)'); set(gca,'FontSize',font_size);
    xticks(xt); xticklabels({'-4','-3','-2','-1','0','1','2','3','4'})
    yticks(zt); 
    yticklabels({'-1.25','-1','-0.75','-0.5','-0.25','-0','0.25',...
        '0.5','0.75','1','1.25','1.5','1.75','2','2.25','2.5','2.75'});
    
    % Visulize the radar heatmap top view
    radar_heatmap_top = squeeze(max(heatmap,[],3));
    subplot(224); imagesc(radar_heatmap_top);
    title("Top View")
    set(gca,'XDir','normal'); set(gca,'YDir','normal');
    colormap jet; caxis([0 maxheat]); colorbar;
    xlabel('x(m)'); ylabel('y(m)'); set(gca,'FontSize',font_size);
    xticks(xt); xticklabels({'-4','-3','-2','-1','0','1','2','3','4'})
    yticks(yt); yticklabels({'3','4','5','6','7','8','9','10'})
end