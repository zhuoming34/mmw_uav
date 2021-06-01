%%% This script is used to read the binary file produced by the DCA1000
%%% and Mmwave Studio
%%% Command to run in Matlab GUI -
close all; clear; clc
%addpath('functions');

Fs = 2047e3;            % Sampling frequency 2047 ksps (kilo samples per second)                    
T = 1/Fs;             % Sampling period       
L = 256;             % Length of signal. Number of samples. 
t = (0:L-1)*T;        % Time 
fc = 60e9;     % why use 60 GHZ? 6843ISK operates at 60-64GHz. 
c = 3e8;
RampEndTime = 133; % us
%sweep_slope = 29.982; % unit MHz/us
sweep_slope = 29.982e+12; % Hz/s

numFrames = 100;
numADCSamples = L;
numChirp = 16; % number of chirps per frame 
numSnap = 64;
numHz = 3; % #horizontal posstions
%SARdata = zeros(numSnap*2,numADCSamples*3); % 256x768, 64*4=256 RX
rf_complex = zeros(numSnap,numADCSamples,4*2*3); % choose only horizontal TXs
%hz_pos = ["1st\", "2nd\"];
%rf_complex_4d = zeros(3,numSnap,numADCSamples,4*2); % temp 4d complex set storing rf raw data from each horizontal position
rf_complex_h1 = zeros(numSnap,numADCSamples,4*2);
rf_complex_h2 = zeros(numSnap,numADCSamples,4*2);
rf_complex_h3 = zeros(numSnap,numADCSamples,4*2);

% find the mean values of corresponding samples from 128 chirps x 64 
% vertical positions x 3 horizontal positions (128*64*3=24576)
svaddr = '..\mmWave\realdata_box\4viewradar\3dinten\';
addr1 = '..\mmWave\realdata_box\4viewradar\';
viewfolder = 'sardata_view'; 
for view = 1:4
    disp(strcat('view: ', num2str(view)))
    for fm = 1:100%numFrames
        disp(strcat('frame: ',num2str(fm)))
        for h = 1:3
            for i = 1:numSnap
                disp(strcat('hz-sn: ', num2str(h),'-',num2str(i)))
                %hz_folder = hz_pos(h);
                hz_folder = h;
                addr2 = i;
                addr3 = '\adc_data.bin';
                addr = strcat(addr1,viewfolder,num2str(view),'\',num2str(hz_folder),'\',num2str(addr2),addr3);
                originData = readDCA1000(addr, numADCSamples); % 4x(32768*3), 32768=128*256

                numSC = numADCSamples*numChirp; % #sample of each frame
                framedata = originData(:,numSC*3*(fm-1)+1:numSC*3*fm); %4x(256*64*3)

                TX1 = framedata(:,numSC*0+1:numSC*1); % 4x256*128
                TX2 = framedata(:,numSC*1+1:numSC*2);
                TX3 = framedata(:,numSC*2+1:numSC*3);
                % TR -> virtual antenna
                % first chirp -> 256x1
                TR11m = reshape(TX1(1,1:numADCSamples),[numADCSamples,1]); 
                TR12m = reshape(TX1(2,1:numADCSamples),[numADCSamples,1]);
                TR13m = reshape(TX1(3,1:numADCSamples),[numADCSamples,1]);
                TR14m = reshape(TX1(4,1:numADCSamples),[numADCSamples,1]);      
                TR31m = reshape(TX3(1,1:numADCSamples),[numADCSamples,1]);
                TR32m = reshape(TX3(2,1:numADCSamples),[numADCSamples,1]);
                TR33m = reshape(TX3(3,1:numADCSamples),[numADCSamples,1]);
                TR34m = reshape(TX3(4,1:numADCSamples),[numADCSamples,1]);
                %{
                TR21m = mean(reshape(TX2(1,:),[numADCSamples,numChirp]),2);
                TR22m = mean(reshape(TX2(2,:),[numADCSamples,numChirp]),2);
                TR23m = mean(reshape(TX2(3,:),[numADCSamples,numChirp]),2);
                TR24m = mean(reshape(TX2(4,:),[numADCSamples,numChirp]),2);

                %}
                % range-fft(sin(phase angle(rawdata))), 256x1
                TR11 = fft(sin(angle(TR11m)));
                TR12 = fft(sin(angle(TR12m)));
                TR13 = fft(sin(angle(TR13m)));
                TR14 = fft(sin(angle(TR14m)));
                TR31 = fft(sin(angle(TR31m)));
                TR32 = fft(sin(angle(TR32m)));
                TR33 = fft(sin(angle(TR33m)));
                TR34 = fft(sin(angle(TR34m)));

                %{ 
                TR21 = fft(sin(angle(TR21m)));
                TR22 = fft(sin(angle(TR22m)));
                TR23 = fft(sin(angle(TR23m)));
                TR24 = fft(sin(angle(TR24m)));
                %}
                %{ 
                %TX1avg = [TR11';TR12';TR13';TR14']; %4x256
                %TX3avg = [TR31';TR32';TR33';TR34'];
                %TRs = [TX1avg' TX3avg']; %256x8, same as below
                %TRs = [TR11 TR12 TR13 TR14 TR21 TR22 TR23 TR24];
                %}
                TRs = [TR11 TR12 TR13 TR14 TR31 TR32 TR33 TR34]; % 256x8


                %rf_complex_4d(h,i,:,:) = TRs;
                switch h
                    case 1
                        rf_complex_h1(i,:,:) = TRs;
                    case 2
                        rf_complex_h2(i,:,:) = TRs;
                    case 3
                        rf_complex_h3(i,:,:) = TRs;
                end
            end
        end

        % flip over elements in each column
        % i-th snapshot is at the bottom
        rf_complex_h1 = flip(rf_complex_h1); % 64x256x8
        rf_complex_h2 = flip(rf_complex_h2); 
        rf_complex_h3 = flip(rf_complex_h3); 

        rf_complex = [rf_complex_h1,rf_complex_h2,rf_complex_h3]; % 64x256x24

        disp('done intiallizing')

        % generate 3d intensity maps
        disp('generate 3d intensity maps')
        radar_heatmap = zeros(256,64,64); % Hawkeye data structure
        %rf_data_full = rf_complex; % full scale
        rf_data_2ss = rf_complex(31:32,:,:); % 2 snapshots
        NTX = size(rf_complex, 3); % 24
        fft_data = zeros(numSnap,numSnap);
        
        for range_idx = 1:numADCSamples
            % fft_data(:,1:NTX) = rf_data_full(:,range_idx,:); % full-scale
            fft_data(1:2,1:NTX) = rf_data_2ss(:,range_idx,:);
            all_fft = fft2(fft_data); % angle fft
            all_fft = fftshift(all_fft);
            %all_fft = flip(flip(all_fft,2));
            temp = abs(all_fft);
            temp2 = temp.'; %(elv,azi) -> (azi,elv)
            radar_heatmap(range_idx,:,:) = temp2;
        end

        save(strcat(svaddr,num2str((fm-1)*4+view),'.mat'), 'radar_heatmap');
        disp(strcat('finished frame: ',num2str(fm))) 
    end
    disp(strcat('finished view: ',num2str(view)))
end

% swra581b.pdf -- Mmwave Radar Device ADC Raw Data Capture
% spruij4a.pdf -- DCA1000EVM Data Capture Card

function [retVal] = readDCA1000(fileName, numADCSamples)
    %% global variables
    % change based on sensor config
    %numADCSamples = 256; % number of ADC samples per chirp
    numADCBits = 16; % number of ADC bits per sample
    numRX = 4; % number of receivers
    numLanes = 2; % do not change. number of lanes is always 2
    isReal = 0; % set to 1 if real only data, 0 if complex data0
    %% read file
    % read .bin file
    fid = fopen(fileName,'r');
    adcData = fread(fid, 'int16');
    % if 12 or 14 bits ADC per sample compensate for sign extension
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
    end
    fclose(fid);
    fileSize = size(adcData, 1);
    %disp(fileSize)
    %disp(size(adcData))
    
    % real data reshape, filesize = numADCSamples*numChirps
    if isReal
        numChirps = fileSize/numADCSamples/numRX;
        LVDS = zeros(1, fileSize);
        %create column for each chirp
        LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
        %each row is data from one chirp
        LVDS = LVDS.';
    else
        % for complex data
        % filesize = 2 * numADCSamples*numChirps
        numChirps = fileSize/2/numADCSamples/numRX;
        %disp("numChirps ="+numChirps)

        LVDS = zeros(1, fileSize/2);
        % Each LVDS lane contains 2 bytes, either 2I or 2Q. 
        %combine real and imaginary part into complex data
        %read in file: 2I is followed by 2Q
        counter = 1;
        for i=1:4:fileSize-1
            LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
            LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
            counter = counter + 2;
        end
        
        %disp("LVDS size ="+size(LVDS));
        
        % create column for each chirp
        LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
        %disp("LVDS="+size(LVDS));
        %each row is data from one chirp
        LVDS = LVDS.';
        %disp("LVDS="+size(LVDS));
    end
    %organize data per RX
    adcData = zeros(numRX,numChirps*numADCSamples);
    for row = 1:numRX
        for i = 1: numChirps
            adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = ...
                LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
        end
    end
    % return receiver data
    %disp("adcData size =" + size(adcData))
    retVal = adcData;
end
