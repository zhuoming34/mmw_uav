% plot point cloud in cartesian coordinate system
function show_ct_pt(ct_coord)
    font_size = 8;
    figure(); scatter3(ct_coord(:,1),ct_coord(:,2),ct_coord(:,3),0.3,'filled','k');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); axis equal;
    xlim([-5 5]); ylim([0 10]); zlim([-1.25,3]); set(gca,'FontSize',font_size);      
    %view(15,30)
end