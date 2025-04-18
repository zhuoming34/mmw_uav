close all; clear; clc;

boxaddr = '../dataset/zmBoxes/';
saveaddr = '../dataset/zmBoxes/';

% ---- load boxes ----
boxname1 = 'box38x59x48'; 
boxname2 = 'box36x58x47'; 
boxname3 = 'box31x50x36'; 
%boxname4 = 'box49x38x58';
box1 = load(strcat(boxaddr,boxname1,'.mat')).myBox; %
box2 = load(strcat(boxaddr,boxname2,'.mat')).myBox; %
box3 = load(strcat(boxaddr,boxname3,'.mat')).myBox; %
%box4 = load(strcat(boxaddr,boxname3,'.mat')).myBox; %

%{
boxt = box2;
figure;
scatter3(boxt(:,1),boxt(:,2),boxt(:,3),0.5,'filled','k');
axis equal;
W = max(boxt(:,1))-min(boxt(:,1));
L = max(boxt(:,2))-min(boxt(:,2));
H = max(boxt(:,3))-min(boxt(:,3));
xlim([-W,W]),ylim([-L,L]),zlim([-H/2,3*H/2]);
%}

% ---- Rotate ---- 
% inline function for 2D rotation
rotate2d =  @(x, M) (x(:, 1:2) * M);
box.rotate = [0,0,0]; % angle in degree, 0~360, clockwise, rotate about origin
%box.rotate = mod(box.rotate,180);
rotate_angle_rad = box.rotate/180*pi;
rotation_matrix1 = [cos(rotate_angle_rad(1)), -sin(rotate_angle_rad(1)); sin(rotate_angle_rad(1)), cos(rotate_angle_rad(1))]; % create rotation matrix
box1(:,1:2) = rotate2d(box1, rotation_matrix1); % rotate the point cloud 

rotation_matrix2 = [cos(rotate_angle_rad(2)), -sin(rotate_angle_rad(2)); sin(rotate_angle_rad(2)), cos(rotate_angle_rad(2))]; % create rotation matrix
box2(:,1:2) = rotate2d(box2, rotation_matrix2); % rotate the point cloud 

rotation_matrix3 = [cos(rotate_angle_rad(3)), -sin(rotate_angle_rad(3)); sin(rotate_angle_rad(3)), cos(rotate_angle_rad(3))]; % create rotation matrix
box3(:,1:2) = rotate2d(box3, rotation_matrix3); % rotate the point cloud 

%{
figure;
scatter3(box1(:,1),box1(:,2),box1(:,3),0.5,'filled','k');
axis equal;
W = max(box1(:,1))-min(box1(:,1));
L = max(box1(:,2))-min(box1(:,2));
H = max(box1(:,3))-min(box1(:,3));
xlim([-W,W]),ylim([-L,L]),zlim([-H/2,3*H/2]);
 %}

% ---- Translation ----
translate_x = -190; translate_y = 0; translate_z = 0; % mm
translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
box1 = box1 + translate; % translate the point cloud

translate_x = 180; translate_y = 0; translate_z = 0; % mm
translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
box2 = box2 + translate; % translate the point cloud

translate_x = -190; translate_y = -45; translate_z = 480; % mm
translate = [translate_x, translate_y, translate_z]; % store translation information in the pc structure
box3 = box3 + translate; % translate the point cloud

% convert unit from mm to m
%box1m = box1/1000; 
%box2m = box2/1000; 
%box3m = box3/1000; 
%box4m = box4/1000; 

% ---- random sampling ----

boxes = [box1;box2;box3]; % mm be saved
%boxes_m = [box1m;box2m;box3m;box4m]; % m be plotted


filename = strcat(saveaddr,'lab_scene2.mat');
save(filename, 'boxes')
disp("saved");

%%
figure;
%scatter3(boxes_m(:,1),boxes_m(:,2),boxes_m(:,3),0.5,'filled','k');
scatter3(boxes(:,1),boxes(:,2),boxes(:,3),0.5,'filled','k');
hold on; scatter3(0,0,1,0.5,'filled','m'); hold off;
axis equal;
%xlim([-4,4]),ylim([-2,8]),zlim([0,3]);
xlim([-4000,4000]),ylim([-2000,8000]),zlim([0,3000]);
xlabel('x (m)'),ylabel('y (m)'),zlabel('z (m)')
            
%%
%{
% ----sensors----
h = 1;
% ---radar---
rhfov_deg = 120; rvfov_deg = 30;
rhfov = rhfov_deg/180*pi; rvfov = rvfov_deg/180*pi;
l1 = h/tan(rvfov/2); l2= 8;
h1 = 0; h2 = l2*tan(rvfov/2);
w = 3; lr = w/tan(rhfov/2);
vfovl1 = [0,0,h;0,l1,h1]; vfovl2 = [0,0,h;0,l2,h2+h];
hfovl1 = [0,0,h;-w,lr,h]; hfovl2 = [0,0,h;w,lr,h];

hres = 2/24; vres = 2/64;
hres_deg = hres/pi*180; vres_deg = vres/pi*180;
hres_p = rhfov/64; vres_p = rvfov/64;
hres_p_deg = rhfov_deg/64; vres_p_deg = rvfov_deg/64;
Nhcell = hres_deg/hres_p_deg; Nvcell = vres_deg/vres_p_deg;

%%
% ---camera---
chfov_deg = 85; cvfov_deg = 54;
chfov = chfov_deg/180*pi; cvfov = cvfov_deg/180*pi;
l3 = h/tan(cvfov/2); l4= 8;
h3 = 0; h4 = l2*tan(cvfov/2);
lc = w/tan(chfov/2);
vfovl3 = [0,0,h;0,l3,h3]; vfovl4 = [0,0,h;0,l4,h4+h];
hfovl3 = [0,0,h;-w,lc,h]; hfovl4 = [0,0,h;w,lc,h];

% ----objs----
blen = 0.5;
hblen = blen/2;
r = 4; % meter

phi_diff = hres*7; 
w_obj = r*sin(phi_diff/2);
l_obj = sqrt((r^2-w_obj^2)-(h-hblen)^2);
theta_diff = asin((h-hblen)/r);%theta_diff = vres*5;
%h_obj = r*sin(theta_diff/2);
h_obj = h-hblen; % hblen*2 + 0.5;
%l_obj2 = sqrt(r^2-h_obj^2);
bsize = [blen blen blen]; % box edge length
obj_srf_ctr = [-w_obj l_obj hblen; w_obj l_obj hblen; -w_obj l_obj h_obj]; % reflecting points
%obj_ctr0 = [-w_obj l_obj+hblen hblen; w_obj l_obj+hblen hblen; -w_obj l_obj+hblen h+h_obj]; 
obj_ctr = obj_srf_ctr + [0 hblen 0; 0 hblen 0; 0 hblen h]; % box's center
obj_org = obj_ctr - bsize/2; % box's origin

% ----plots----
figure(1); 
p1 = scatter3(0,0,h,50,'filled','m'); % sensor = [0,0,1];

hold on; p2 = plot3(hfovl1(:,1)',hfovl1(:,2)',hfovl1(:,3)','-.b'); % radar FOV
hold on; plot3(hfovl2(:,1)',hfovl2(:,2)',hfovl2(:,3)','-.b'); 
hold on; p3 = plot3(vfovl1(:,1)',vfovl1(:,2)',vfovl1(:,3)','--b');
hold on; plot3(vfovl2(:,1)',vfovl2(:,2)',vfovl2(:,3)','--b'); 

hold on; p4 = plot3(hfovl3(:,1)',hfovl3(:,2)',hfovl3(:,3)','-.r'); % camera FOV
hold on; plot3(hfovl4(:,1)',hfovl4(:,2)',hfovl4(:,3)','-.r'); 
hold on; p5 = plot3(vfovl3(:,1)',vfovl3(:,2)',vfovl3(:,3)','--r');
hold on; plot3(vfovl4(:,1)',vfovl4(:,2)',vfovl4(:,3)','--r'); 

plotcube(bsize, obj_org(1,:), 0.5, [1 0 0]); % object 1
plotcube(bsize, obj_org(2,:), 0.5, [0 1 0]); % object 2
plotcube(bsize, obj_org(3,:), 0.5, [0 0 1]); % object 3
plotcube(bsize, [0-0.25 5.5 1-0.25], 0.5, [0 1 1]); % object 4

plotcube([10 10 0.001], [-3 -3 -0.001], 0.1, [0 0 0]); % floor

legend([p1 p2 p3 p4 p5],{'Sensors','Radar hFOV','Radar vFOV','Camera hFOV','Camera vFOV'})
axis equal;
xlabel('x(m)'),ylabel('y(m)'),zlabel('z(m)')
xlim([-2.5,2.5]),ylim([-2.5,6.5]),zlim([0,3])

function plotcube(varargin)
% PLOTCUBE - Display a 3D-cube in the current axes
%
%   PLOTCUBE(EDGES,ORIGIN,ALPHA,COLOR) displays a 3D-cube in the current axes
%   with the following properties:
%   * EDGES : 3-elements vector that defines the length of cube edges
%   * ORIGIN: 3-elements vector that defines the start point of the cube
%   * ALPHA : scalar that defines the transparency of the cube faces (from 0
%             to 1)
%   * COLOR : 3-elements vector that defines the faces color of the cube
%
% Example:
%   >> plotcube([5 5 5],[ 2  2  2],.8,[1 0 0]);
%   >> plotcube([5 5 5],[10 10 10],.8,[0 1 0]);
%   >> plotcube([5 5 5],[20 20 20],.8,[0 0 1]);
% Default input arguments
inArgs = { ...
  [10 56 100] , ... % Default edge sizes (x,y and z)
  [10 10  10] , ... % Default coordinates of the origin point of the cube
  .7          , ... % Default alpha value for the cube's faces
  [1 0 0]       ... % Default Color for the cube
  };
% Replace default input arguments by input values
inArgs(1:nargin) = varargin;
% Create all variables
[edges,origin,alpha,clr] = deal(inArgs{:});
XYZ = { ...
  [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
  [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
  [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
  [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
  [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
  [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
  };
XYZ = mat2cell(...
  cellfun( @(x,y,z) x*y+z , ...
    XYZ , ...
    repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
    repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
    'UniformOutput',false), ...
  6,[1 1 1]);
cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
  repmat({clr},6,1),...
  repmat({'FaceAlpha'},6,1),...
  repmat({alpha},6,1)...
  );
view(3)
end
%}