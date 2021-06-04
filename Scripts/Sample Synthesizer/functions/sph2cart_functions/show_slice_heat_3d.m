% plot 3d heatmap slice
function show_slice_heat_3d(heatmap_ct)
    m = size(heatmap_ct, 2); % az
    n = size(heatmap_ct, 3); % el
    l = size(heatmap_ct, 1); % rg
    xi = linspace(1,l,l); % range
    yi = linspace(1,m,m); % azimuth
    zi = linspace(1,n,n); % elevation
    [XX,YY,ZZ] = meshgrid(yi,xi,zi); % [l * m * n]

    xslice = 1:m;    % location of y-z planes
    yslice = 1:l;    % location of x-z plane
    zslice = 1:n;    % location of x-y planes
    
    figure(); h = slice(XX,YY,ZZ,heatmap_ct,xslice,yslice,zslice);
    xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
    xlim([0 64]); ylim([0 256]); zlim([0,64]);
    xt = linspace(0,64,9); xticks(xt); 
    xticklabels({'-4','-3','-2','-1','0','1','2','3','4'})
    yt = linspace(0,256,8); yticks(yt); 
    yticklabels({'3','4','5','6','7','8','9','10'})
    zt = linspace(0,64,17); zticks(zt); 
    zticklabels({'-1.25','-1','-0.75','-0.5','-0.25','-0','0.25',...
        '0.5','0.75','1','1.25','1.5','1.75','2','2.25','2.5','2.75'});
    %set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp');
    set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','flat');
    %set(h,'EdgeColor','none','FaceColor','flat','FaceAlpha','0.01');
    % set transparency to correlate to the data values.
    alpha('color'); %alpha(h, 0.01);
    colorbar; colormap jet;
end