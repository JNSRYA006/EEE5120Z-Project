radar = radarParameters.radar1;

first_intercept = radar.TOA;
AOA = radar.AOA;
PRI = radar.PRI;
beamwidth = radar.Beamwidth; 
scanPeriod = radar.scanPeriod; % s
scanRate = 360/scanPeriod;
fprintf('%.0f deg/s\n',scanRate);
sensitivity = -50; % dBm
%% Testing List generation methods
observationWindow = 3; % s
detection = false;
num_pulses = floor((observationWindow-first_intercept)/PRI);
fprintf('Number of pulses of observation window of %.0fs is: %.0f\n',observationWindow,num_pulses);

for i = 1:num_pulses
    TOA = (i-1)*PRI+first_intercept;
    TOA_dig = discretiseVal(TOA,PDWParameters.TOA.min,PDWParameters.TOA.resolution,PDWParameters.TOA.max-PDWParameters.TOA.min);
    scanRate = 360/scanPeriod;
    AOA_offset = scanRate*TOA;
    angle = mod(radar.AOA+AOA_offset,360);
    disp(angle);

    amplitude = abs(sind((angle*pi)/2*beamwidth));
    amplitude_db = 20*log10(amplitude) + radar.peakAmp;

    if amplitude_db >= sensitivity
        if ~detection
            detection_angle_first = angle;
            detection = true;
            detection_angle_last = angle + 2*beamwidth; %wrong. Shouldn't be that wide
        end
    else
        amplitude_dB = -50;
    end

    if (angle >= detection_angle_first && angle <= detection_angle_last && amplitude_db >= sensitivity)
        amplitude_dB = amplitude_db;
    else
        amplitude_dB = -50;
    end
    PDW(i) = generatePDW(PDWParameters,TOA,angle,amplitude_dB,radar.Freq,radar.PW,radarParameters);
    % PDW.
    disp(PDW);
end    

%% Plotting
AOA_plot = zeros(1,num_pulses);
Amp_plot = zeros(1,num_pulses);
for i = 1:num_pulses
    AOA_plot(i) = PDW(i).AOA;
    Amp_plot(i) = PDW(i).Amp;
end
figure;
plot(AOA_plot,Amp_plot,'-o');
xlabel('Angle of Arrival [deg]');
ylabel('Amplitude [dB]');
ylim([-50 -30]);
xlim([AOA-3*beamwidth AOA+3*beamwidth]);
title('Amplitude vs. AOA');
grid on;
%% Playing around

n = 1000;
% Define the angle theta from 0 to 360 degrees
theta = linspace(0, 360, n);
theta_mainlobe = theta;
theta_mainlobe(theta_mainlobe>2*beamwidth) = 0;

% Compute the function f(theta)
amplitude = sind((theta*pi)/(2*beamwidth));

% Compute the amplitude in dB
amplitude_dB = 20*log10(amplitude)+(radar.peakAmp);

% Num of pulses per beamwidth rotation
t_beamwidth_rotation = beamwidth/scanRate; % Time for 3dB beamwidth to rotate through
t_mainlobe_rotation = 2*beamwidth/scanRate; % Time for mainlobe to rotate through full size
fprintf('%.3f ms for rotation through %.2f deg (mainlobe)\n',t_mainlobe_rotation*1e3,2*beamwidth);
PRI = radar.PRI;
pulse_beamwidth_rotation = floor(t_beamwidth_rotation/PRI); % Number of pulses emitted for 3dB beamwidth to rotate
pulse_mainlobe_rotation = floor(t_mainlobe_rotation/PRI); % Number of pulses emitted for whole mainlobe to rotate
fprintf('%.0f full pulses emitted in rotation through beamwidth\n',pulse_mainlobe_rotation);
disp('-----------------------------------------------');

boresight_time = (beamwidth/360)*scanPeriod;

t_first = scanTimeBeamwidth(detection_index_first);
% disp(t_first*1e3);
t_detect = t_mainlobe_rotation - 2*first_intercept; % Time for which mainlobe is detected (> sensitivity of ESM receiver)
fprintf('Radar detected for %.5f ms\n',t_detect*1e3);
mainlobe_detected_pulses = floor(t_detect/PRI); % Pulses detected in mainlobe based on time it is detected for (not accounting for any offsets)
fprintf('This means that %.0f pulses are detected in the main lobe\n',mainlobe_detected_pulses);
% disp(t_detect*1e3);
% first_intercept = radarParameters.radar3.TOA; % Actually t_first based on Table 2, because that is first sensitivity crossing
last_intercept = first_intercept + t_detect;
fprintf('First intercept: %.4f; Last intercept: %.4f (ms)\n',first_intercept*1e3,last_intercept*1e3);
fprintf('Over the entire mainlobe, %.0f pulses are transmitted. ',pulse_mainlobe_rotation);
fprintf('Of these, %.0f, are detected based on t_detect.\n',floor((t_mainlobe_rotation-PRI)/PRI));
disp('-----------------------------------------------');

fprintf('Initial angle (AOA Tab. 2) is: %.0f at t = %0.4f us\n',radar.AOA,radar.TOA*1e3);
disp('This is the first detection (i.e. amplitude_db >= sensitivity)');
fprintf('After this detection, ')


%% Functions
function discretised = discretiseVal(value,min_val,resolution,range)
    norm_val = (value-min_val)/range;
    levels = range/resolution;
    discrete_level = round(norm_val*levels);
    discretised = min_val + discrete_level*resolution;
    % range = (2^bits -1)*resolution;
end


