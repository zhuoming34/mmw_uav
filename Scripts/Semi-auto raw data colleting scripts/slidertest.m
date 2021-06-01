%% test Arduino Uno and Adafruit motor shield V2 with stepper motors
close all; clear; clc;
a = arduino('COM15','Uno','Libraries','Adafruit\MotorShieldV2');
shield = addon(a,'Adafruit\MotorShieldV2');
%addrs = scanI2CBus(a,0);
% setup stepper motors
sm1 = stepper(shield,1,200,'RPM',100,'StepType','double'); % slider 1, M1M2
sm2 = stepper(shield,2,200,'RPM',500,'StepType','double'); % slider 2, M3M4
disp(a),disp(shield),disp(sm1),disp(sm2);%disp(addrs);

%% 
move(sm2,31500); % 500 = 2.5mm
disp('moved')
release(sm2);
disp('released')
%%
move(sm1,-8000);
release(sm1);

%pause(1); % in second



