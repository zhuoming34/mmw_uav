% scale the real depth images to a specified boundary and size
close all; clear; clc;
addr = '..\mmWave\realbox0526\20210526\cam';
svaddr = '..\mmWave\realbox0526\';


numc = 1280; numr = 720;
r1 = 0.1*1000; r2 = 15.0*1000; % meter, min/max range
%N_col_l = 407; %N_col_r = 407; %N_row_t = 244; %N_row_b = 127; 

f_px = 700; % focal lenght in pixel
pxSize = 0.004; % mm
f = f_px*pxSize/1000; % 2.8mm -> 0.0028m
% projection plane 
ppH = numc*pxSize/1000; ppV = numr*pxSize/1000; % in m
xl = -ppH/2; xr = ppH/2; zb = -ppV/2; zt = ppV/2;

scene_corners = [-1,3,-0.5; -1,3,1;  1,3,-0.5; 1,3,1]; % the space selected for scaling depth images
sn_x = scene_corners(:,1); sn_y = scene_corners(:,2); sn_z = scene_corners(:,3);
% scene projection
sn_py = linspace(f,f,size(sn_y,1))'; % focal length or y values of projefction plane
sn_ratio = sn_y./sn_py; 
sn_px = -sn_x./sn_ratio;
sn_pz = -sn_z./sn_ratio;

N_col_l = ceil(abs((sn_px(3)-xl))/(pxSize/1000));
N_col_r = ceil(abs((sn_px(1)-xr))/(pxSize/1000));
N_row_t = ceil(abs((sn_pz(1)-zt))/(pxSize/1000));
N_row_b = ceil(abs((sn_pz(2)-zb))/(pxSize/1000));

map1 = jet; map2 = gray;

%%
for idx = 0:99
    for cam = 1:4
        addr2 = strcat(addr,num2str(cam),'\');
        
        if idx < 10
            matname_prefix = 'depthmat00000';
            imgname_prefix = 'depthimg00000';
        else
            matname_prefix = 'depthmat0000';
            imgname_prefix = 'depthimg0000';
        end
        filename_mat = strcat(addr2,matname_prefix,num2str(idx),'.mat');
        depthmat = load(filename_mat).depmat;
        %filename_img = strcat(addr2,'left00000',num2str(idx),'.png');
        %orgImg = imread(filename_img);
        filename_img2 = strcat(addr2,imgname_prefix,num2str(idx),'.png');
        orgImg2 = imread(filename_img2);
               
        %depthmat(depthmat <= r1) = r1;
        %depthmat(depthmat >= r2) = r2;
        
        depthmat(depthmat <= r1) = r1;
        %depthmat(depthmat >= r2) = r2;
        depthmat(1:295,:) = r2;
        depthmat(460:end,:) = r2;
        depthmat(:,1:550) = r2;
        depthmat(:,780:end) = r2;
             
        switch cam
            case 1
                depthmat(1:360,1:650) = r2;
                depthmat(:,1:565) = r2;
                depthmat(:,765:end) = r2;
                %depthmat(depthmat >= 3.5e+03) = r2;
            case 2
                depthmat(1:360,630:end) = r2;
                depthmat(:,675:end) = r2;
                %depthmat(depthmat >= 4e+03) = r2;
            case 3
                depthmat(1:370,675:end) = r2;
                depthmat(:,1:585) = r2;
                depthmat(440:end,:) = r2;
                depthmat(:,750:end) = r2;
                %depthmat(depthmat >= 3.5e+03) = r2;
            case 4
                depthmat(1:360,1:700) = r2;
                depthmat(1:365,765:end) = r2;
                depthmat(:,1:655) = r2;
                depthmat(1:300,:) = r2;
                %depthmat(depthmat >= 3e+03) = r2;
        end
        
        depthmat(1,numc/2) = r1;
        
        % normalize pixel values to [0,1]
        I0 = depthmat;
        I = ( I0 - min(min(I0)) ) ./ ( max(max(I0)) - min(min(I0)) );
        
        % extend and crop to the same size of scene  
        [Ib,N_row,N_col] = depScale(I,sn_px,xl,xr,zt,zb);
        I2 = Ib;
        
        % saving images
        Dep1 = abs(I - 1)*255;
        Dep2 = abs(I2 - 1)*255;
        reszImg = imresize(Dep2, [128,128]);

        % idx starts from 0
        newidx = idx+1;%*4+cam + 2400;
        %outputFileName0 = strcat(svaddr,'original2/',num2str(newidx),'.png');
        %imwrite(Dep1, map2, outputFileName0); 

        outputFileName1 = strcat(svaddr,'original/cam',num2str(cam),'/',num2str(newidx),'.png');
        imwrite(Dep1, map2, outputFileName1); 

        outputFileName2 = strcat(svaddr,'extend/cam',num2str(cam),'/',num2str(newidx),'.png');
        imwrite(Dep2, map1, outputFileName2); 

        outputFileName3 = strcat(svaddr,'color128/cam',num2str(cam),'/',num2str(newidx),'.png');
        imwrite(reszImg, map1, outputFileName3);

        outputFileName4 = strcat(svaddr,'gray128/cam',num2str(cam),'/',num2str(newidx),'.png');
        imwrite(reszImg, map2, outputFileName4); 
        
        disp(strcat(num2str(idx)," - ", num2str(cam)));
    end
end
disp('finished')
%%
show_photo = 0;
if show_photo == 1
    %imshow(orgImg, map);
    %figure(1)
    %heatmap(orgImg(:,:,1));
    %colormap jet;
    %grid off;
    figure
    heatmap(orgImg2(:,:,1));
    colormap jet;
    grid off;
end

function [Ib,N_row,N_col] = depScale(I,sn_px,xl,xr,zt,zb)        
    if sn_px(3) <= xl || sn_px(1) >= xr
        col_l = ones(numr,N_col_l); col_r = ones(numr,N_col_r); % extend
        Ic = [col_l,I,col_r];
        N_col = numc + N_col_l + N_col_r;
    else
        Ic = I(:,N_col_l+1:end-N_col_r); % crop
        N_col = numc - N_col_l - N_col_r;
    end

    if sn_pz(1) >= zt
        row_t = ones(N_row_t,N_col); % extend
        It = [row_t;Ic];
        N_row = numr + N_row_t;
    else
        It = Ic(N_row_t+1:end,:); % crop
        N_row = numr - N_row_t;
    end

    if sn_pz(2) <= zb
        row_b = ones(N_row_b,N_col); % extend
        Ib = [Ic;row_b];
        N_row = N_row + N_row_b;
    else
        Ib = It(1:end-N_row_b,:); % crop
        N_row = N_row - N_row_b;
    end
end