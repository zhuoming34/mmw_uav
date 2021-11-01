function [depImg,reszImg,squrImg] = pc2dep(ptCloud)
% converts a point cloud into a 2d depth image
% Input: 'ptCloud': a point cloud of a surface
% Output: 'depImg': original size of a depth image
%         'reszImg': extended/cropped depth image
%         'squrImg': scaled reszImg to be square

    variable_library_camera; % load camera configurations
    
    x = ptCloud(:,1); y = ptCloud(:,2); z = ptCloud(:,3);
    
    % perspective projection for each point
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
    xl = -ppH/2; xr = ppH/2; zb = -ppV/2; zt = ppV/2;
    xx = linspace(xr,xl,numc); zz = linspace(zb,zt,numr);
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
    % usually N_col_l = N_col_r
    N_col_l = ceil(abs((sn_px(3)-xl))/(pxSize/1000));
    N_col_r = ceil(abs((sn_px(1)-xr))/(pxSize/1000));
    N_row_b = ceil(abs((sn_pz(1)-zt))/(pxSize/1000));
    N_row_t = ceil(abs((sn_pz(2)-zb))/(pxSize/1000));
    if sn_px(3) <= xl || sn_px(1) >= xr
        col_l = ones(numr,N_col_l); col_r = ones(numr,N_col_r); % extend
        Ic = [col_l,I,col_r];
        N_col = numc + N_col_l + N_col_r;
    else
        Ic = I(:,N_col_l+1:end-N_col_r); % crop
        N_col = numc - N_col_l - N_col_r;
    end

    if sn_pz(1) >= zt % positive on projection plane
        row_b = ones(N_row_b,N_col); % extend
        Ib = [Ic;row_b];
        N_row = numr + N_row_b;
    else
        Ib = Ic(1:end-N_row_b,:); % crop
        N_row = numr - N_row_b;
    end

    if sn_pz(2) <= zb % negative on projection plane
        row_t = ones(N_row_t,N_col); % extend
        It = [row_t;Ic];
        N_row = N_row + N_row_t;
    else
        It = Ib(N_row_t+1:end,:); % crop
        N_row = N_row - N_row_t;
    end
    I2 = It;

    % saving images
    depImg = abs(I - 1)*255;
    reszImg = abs(I2 - 1)*255;
    squrImg = imresize(reszImg, [128,128]);
    %{
    map = jet;
    map2 = gray;

    outputBaseFileName1 = strcat(outaddr1,'cam',num2str(cam+camos),'/',num2str(idx+os),'.png');
    imwrite(Dep1, map, outputBaseFileName1); 

    outputBaseFileName2 = strcat(outaddr2,'cam',num2str(cam+camos),'/',num2str(idx+os),'.png');
    imwrite(Dep2, map, outputBaseFileName2); 

    outputBaseFileName3 = strcat(outaddr3,'cam',num2str(cam+camos),'/',num2str(idx+os),'.png');
    imwrite(reszImg, map, outputBaseFileName3);

    outputBaseFileName4 = strcat(outaddr4,'cam',num2str(cam+camos),'/',num2str(idx+os),'.png');
    imwrite(reszImg, map2, outputBaseFileName4); 
    %}
end