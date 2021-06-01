%%% SAR automation ver.4
%%% after configuration in mmWave Studio, automatically click on bottons to
%%% arm the DCA1000 and trigger frames, move slider, and save files.
close all; clear; clc;

%% IMPORTANT NOTE

% Check your COM port# for Arduino

% !!!
% Manually trigger one frame first (DCA arm -> Trig frame -> PostProc)
% Make sure the window of post-processing result is pop-up (takes a while) 
% before running the following data-collecting process,
% so that the post processing will be executed right after the button was pressed
% and the pop-up window will not block the mouse pointer event. 

% Use get(0,'PointerLocation') command to get the locations of buttons you
% are going to press (DCA Arm, Trig Frame, PostProc) for your screen.
% Modify the corresponding values below.

% Modify the time for receving the reflectd signal (see line 103 below)

% For first-time use, run it section by section to see how everything works

%% set up file directory
N_ss = 64; % Number of snapshots
adcFolder = 'C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc';
adcFile = strcat(adcFolder,'\adc_data.bin');
rootFolder = '..\mmWave\SAR automation';
dataFolder = 'sardata';
cd(rootFolder)
% create folders to store SAR data
if ~exist(dataFolder, 'dir')
       mkdir(dataFolder) 
end

%% set up mouse object and event
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
screenSize = get(0, 'screensize');

%% arduino motor driver
a = arduino('COM15','Uno','Libraries','Adafruit\MotorShieldV2');
shield = addon(a,'Adafruit\MotorShieldV2');
%addrs = scanI2CBus(a,0);
% setup stepper motors
sm1 = stepper(shield,1,200,'RPM',100,'StepType','double'); % slider 1, M1M2, LeftRight
sm2 = stepper(shield,2,200,'RPM',500,'StepType','double'); % slider 2, M3M4, UpDown
disp(a),disp(shield),disp(sm1),disp(sm2);%disp(addrs);

%% collecting process
up_dn = -1; % -1 => towards the stepper motor, 1 => away. motor of the vertical slide is on the top 
start_pos = 1; end_pos = N_ss; % at first horizontal posiiton, starts from bottom to top
for hz = 1:3
    disp(strcat('start from: ',num2str(start_pos),' end up to:', num2str(end_pos)))
    for ss = start_pos:up_dn*(-1):end_pos
        indFolder = strcat(rootFolder,'\',dataFolder,'\',num2str(hz),'\',num2str(ss));
        if ~exist(indFolder, 'dir')
           mkdir(indFolder) 
        end

        % trigger radar
        % ErrStatus =RtttNetClientAPI.RtttNetClient.SendCommand(Lua_String);
        % use -> get(0,'PointerLocation') <- command to obtain the coordinates of pixel to click on
        % Arm DCA1000: (240,598) -> (240,screenSize(4)-598)
        %pause(1)
        disp(strcat("Capturing position (",num2str(hz),', ',num2str(ss),')'));
        mouse.mouseMove(240,screenSize(4)-598); % (horizontal,verital)
        pause(0.5)
        mouse.mousePress(InputEvent.BUTTON1_MASK);
        mouse.mouseRelease(InputEvent.BUTTON1_MASK);
        % trigger frame: (330,598) -> (330,screenSize(4)-598)
        pause(1)
        mouse.mouseMove(330,screenSize(4)-598); % (horizontal,verital)
        pause(0.5)
        mouse.mousePress(InputEvent.BUTTON1_MASK);
        mouse.mouseRelease(InputEvent.BUTTON1_MASK);

        pause(3.5) % wait for radar signal: 
        % for a 49.7% duty cycle
        % 128 chirp/frame -> 180ms/frame * 100frames = 18s 
        % 64 chirp/frame -> 90ms/frame * 100frames = 9s
        % 16 chirp/frame -> 22.5ms/frame * 100frames = 2.25s
        
        % post processing: (620,598) -> (330,screenSize(4)-598)
        pause(1)
        mouse.mouseMove(620,screenSize(4)-598); % (horizontal,verital)
        pause(0.5)
        mouse.mousePress(InputEvent.BUTTON1_MASK);
        mouse.mouseRelease(InputEvent.BUTTON1_MASK);
        
        next_pos = ss + up_dn*(-1);
        if (up_dn == -1) && (ss == N_ss) 
            disp(strcat("stay at current vertical position (",num2str(hz),', ',num2str(ss),')'));  
            next_pos = ss;
        elseif (up_dn == 1) && (ss == 1) 
            disp(strcat("stay at current vertical position (",num2str(hz),', ',num2str(ss),')'));
            next_pos = ss;
        else
            % move slider vertically
            disp(strcat("Moving to position (",num2str(hz),', ',num2str(next_pos),')'));
            move(sm2,500*up_dn); % about 2s
            release(sm2);         
            pause(1);
        end
        
        % copy file
        copyfile(adcFile,indFolder)
        % done
        disp(strcat("Saved position (",num2str(hz),', ',num2str(ss),')'));
    end
    if hz < 3
        % move slider horizontally
        disp(strcat("Moving to next horizontal position (",num2str(hz+1),', ',num2str(next_pos),')'));
        move(sm1,4000); % about 16ss
        release(sm1);
        %disp(strcat("Moved to position (",num2str(hz),', ',num2str(ss),')'));
        up_dn = up_dn * (-1); % change vertical direction
        old_start_pos = start_pos; % switching index
        start_pos = end_pos; end_pos = old_start_pos;
        pause(1)
    end    
end
disp("finished")
disp("Reseting slider position")
move(sm2,31500); % 500 = 2.5mm
release(sm2);
move(sm1,-8000);
release(sm1);
disp("done")
