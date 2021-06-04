% construct a set of grid center points 
function [ptGrid,ptGrid_heat] = gridpts_contruct(N_x,N_y,N_z,scene_lim)
    xs = linspace(scene_lim(1,1),scene_lim(1,2),N_x);
    ys = linspace(scene_lim(2,1),scene_lim(2,2),N_y);
    zs = linspace(scene_lim(3,1),scene_lim(3,2),N_z);
    W = length(xs); L = length(ys); H = length(zs);
    numPoint = W*L*H;
    ptGrid = zeros(numPoint,3);
    ptGrid_heat = zeros(numPoint,1);
    ptgrid_idx = 1;
    for i = 1:N_x
        for j = 1:N_y
            for k = 1:N_z
                ptGrid(ptgrid_idx,:) = [xs(i),ys(j),zs(k)];  
                ptgrid_idx = ptgrid_idx + 1;
            end
        end
    end
end