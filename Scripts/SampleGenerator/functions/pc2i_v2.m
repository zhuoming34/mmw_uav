%% modified function
% 01/09/2021
% added range limits of depth camera as input 
function I = pc2i_v2(x,y,z,numr,numc,lmin,lmax)
% By: Vahid Behravan
% This function converts a point cloud (given in x,y,z) to a gray scale image
% We assume the ToF camera is alligned with 'y-axis' 
%
% x,y,z: coordinate vectors of all points in the cloud
% numr: desired number of rows of output image
% numc: desired number of columns of output image
% I   : output gray scale image
%
% Example useage:
%   I = pointcloud2image( x,y,z,250,250,0.1,15);
%   figure;  imshow(I,[]);
%-----------------------------------------------------------
% depth calculation
d = sqrt( x.^2 + y.^2 + z.^2);
% grid construction
xl = min(x); xr = max(x); zl = min(z); zr = max(z);
xx = linspace(xl,xr,numc); zz = linspace(zl,zr,numr);
[X,Z] = meshgrid(xx,zz);
grid_centers = [X(:),Z(:)];
% classification
clss = knnsearch(grid_centers,[x,z]); 
% defintion of local statistic
local_stat = @(x)mean(x);
% data_grouping
class_stat = accumarray(clss,d,[numr*numc 1],local_stat);
% 2D reshaping
class_stat_M  = reshape(class_stat , size(X)); 
% adding range limits at front and back
s1 = size(class_stat_M,1);
s2 = size(class_stat_M,2);
class_stat_M(1,s2/2) = lmin;
class_stat_M(s1,s2/2) = lmax;
% Force un-filled cells to the brightest color
class_stat_M (class_stat_M == 0) = max(max(class_stat_M));
% flip image horizontally and vertically
%I = class_stat_M(end:-1:1,end:-1:1);
I0 = class_stat_M(end:-1:1,:);
% normalize pixel values to [0,1]
I = ( I0 - min(min(I0)) ) ./ ( max(max(I0)) - min(min(I0)) );
end