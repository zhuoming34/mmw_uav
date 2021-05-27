% Single box generator, 2020/11/29
% Zhuoming Huang 2020-21 Senior Design
% Run it section by section to prevent stalling
close all; clear; clc;

saveaddr = '../dataset/zmBoxes/';

W = 310; % mm
L = 500;
H = 360;
step = 1; %mm

%x = (1:step:L)-ceil(L/2);
%y = (1:step:W)-ceil(W/2);
z = (0:step:H);
x = (-W/2:step:W/2);%(0:step:L);
y = (-L/2:step:L/2);%(0:step:W);

numPoint = 2*W*H + 2*(L-2)*H + 2*(L-2)*(W-2);
myBox = zeros(numPoint,3);
counter = 0;

showbox = 0; % flag to control showing plot or not

%% 
% need the surface only
for i = [x(1),x(end)]
    for j = y(1):step:y(end)
        for k = z(1):step:z(end)
            counter = counter + 1;
            myBox(counter,:) = [i,j,k];
            
        end
    end
end

%%
for j = [y(1),y(end)]
    for i = x(2):step:x(end-1)
        for k = z(1):step:z(end)
            counter = counter + 1;
            myBox(counter,:) = [i,j,k];
        end
    end
end
%%
for k = [z(1),z(end)]
    for i = x(2):step:x(end-1)
        for j = y(2):step:y(end-1)
            counter = counter + 1;
            myBox(counter,:) = [i,j,k];
        end
    end
end
disp("finished");

% save points
filename = strcat(saveaddr,'box',num2str(W/10),'x',num2str(L/10),'x',num2str(H/10),'.mat');
save(filename, 'myBox')
disp("saved");

% plot
if showbox == 1
    figure;
    scatter3(myBox(:,1),myBox(:,2),myBox(:,3),0.5,'filled','k');
    xlabel('x(mm)'),ylabel('y(mm)'),zlabel('z(mm)');
    axis equal;
    xlim([-W,W]),ylim([-L,L]),zlim([-H/2,3*H/2]); 
end

