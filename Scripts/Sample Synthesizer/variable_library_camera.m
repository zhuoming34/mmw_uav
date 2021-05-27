% Generate 2d depth image from camera reflector point cloud with perspective projection
% https://support.stereolabs.com/hc/en-us/articles/360007395634-What-is-the-camera-focal-length-and-field-of-view-
% zed mini image output resolution of HD720
% width = 1280, height = 720, pixel number = 1280*720
% has a focal length of 700 pixels
% pixel size of 0.004mm
% field of view: fov = 2*arctan(pixelNumber/(2*focalLength))*(180/pi)
% VFOV = 54deg, HFOV = 85deg
% "The focus distance is fixed. All objects at distances from 28cm out to infinity will be sharp"
% depth range = [0.1, 15] meters

numc = 1280; numr = 720; % num of columns and rows
f_px = 700; % focal lenght in pixel
pxSize = 0.004; % pixle size: mm
f = f_px*pxSize/1000; % focal length in meter: 2.8mm -> 0.0028m

% projection plane size in m
ppH = numc*pxSize/1000; ppV = numr*pxSize/1000; 

% horizontal/vertical field of view
hfov = 2*atan(numc/2/f_px); vfov = 2*atan(numr/2/f_px); % rad
hfov_deg = hfov/pi*180; vfov_deg = vfov/pi*180;

% sensing range
r1 = 0.1; r2 = 15.0; % meter, min/max range/depth
h = r2*tan(hfov/2); v = r2*tan(vfov/2); % max angle the camera covers

% scene boundary
scene_corners = [-4,3,-1.25; -4,3,2.75;  4,3,-1.25; 4,3,2.75];
%pts = [ptCloud;scene_corners];
sn_x = scene_corners(:,1); sn_y = scene_corners(:,2); sn_z = scene_corners(:,3);

% scene projection
sn_py = linspace(f,f,size(sn_y,1))'; % focal length
sn_ratio = sn_y./sn_py;
sn_px = -sn_x./sn_ratio;
sn_pz = -sn_z./sn_ratio;

