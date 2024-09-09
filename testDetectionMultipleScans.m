%% Initialise parameters
radar = radarParameters.radar1;

first_intercept = radar.TOA;
AOA = radar.AOA;
PRI = radar.PRI;
beamwidth = radar.Beamwidth; 
mainlobe = 2*beamwidth;
scanPeriod = radar.scanPeriod; % s
scanRate = 360/scanPeriod;
sensitivity = -50; % dBm
observationWindow = 3; % s

%% Calculate Tx Mainlobe

% Num of pulses per beamwidth rotation
disp('-----------------------------------------------');
fprintf('Scan rate of Tx radar: %.0f deg/s\n',scanRate);
t_beamwidth_rotation = beamwidth/scanRate; % Time for 3dB beamwidth to rotate through
t_mainlobe_rotation = mainlobe/scanRate; % Time for mainlobe to rotate through full size
fprintf('%.3f ms for rotation through %.2f deg (mainlobe)\n',t_mainlobe_rotation*1e3,mainlobe);
pulse_beamwidth_rotation = floor(t_beamwidth_rotation/PRI); % Number of pulses emitted for 3dB beamwidth to rotate
pulse_mainlobe_rotation = floor(t_mainlobe_rotation/PRI); % Number of pulses emitted for whole mainlobe to rotate
fprintf('%.0f full pulses emitted in rotation through full mainlobe\n',pulse_mainlobe_rotation);
disp('-----------------------------------------------');

n = 100000;
% Define the angle theta from 0 to 360 degrees
theta = linspace(0, 360, n);
% theta_mainlobe = theta;
% theta_mainlobe(theta_mainlobe>mainlobe) = 0;
theta_mainlobe = linspace(0,mainlobe,n);

% Compute the function f(theta)
amplitude = sin((theta_mainlobe*pi)/(2*beamwidth));

% Compute the amplitude in dB
amplitude_dB = 20*log10(amplitude)+(radar.peakAmp);

% Find indices where amplitude is above sensitivity of ESM receiver
above_sensitivity = find(amplitude_dB > sensitivity);

% Update amplitude and theta arrays to match these indices
amplitude_dB_above_sensitivity = amplitude_dB(above_sensitivity);
theta_mainlobe_above_sensitivity = theta_mainlobe(above_sensitivity);

% Find the angle at which the sensitivity is met
theta_first_detection = theta_mainlobe_above_sensitivity(1);
amplitude_first_detection = amplitude_dB_above_sensitivity(1);
fprintf('First detection occurs at: %.4f degrees with an amplitude of %.4f dB.\n', ...
    theta_first_detection,amplitude_first_detection); % This will vary 
% slightly based on the resolution of the theta linspace

% Need to offset this first detection to occur at t = TOA (Based on Table
% 2)
% Fist detect per Table 2
esm_first_detection = radar.TOA;

% All scan time
scan_time = linspace(0,scanPeriod,n);
% Mainlobe scan time
scan_time_mainlobe = linspace(0,t_mainlobe_rotation,n);
t_boresight = (beamwidth/360)*scanPeriod;

% Calculate shift so first intercept happens at esm_first_detection
time_first_detection = scan_time_mainlobe(above_sensitivity(1));
time_shift = esm_first_detection - time_first_detection;

% Shift scan_time array
scan_time_shift = scan_time+time_shift;
scan_time_mainlobe_shift = scan_time_mainlobe+time_shift;

% Calculate the number of pulses above the sensitivity threshold
t_pulses_detected = scan_time_mainlobe(above_sensitivity(end))-scan_time_mainlobe(above_sensitivity(1));
pulse_detection_rotation = floor(t_pulses_detected/PRI); % Number of pulses emitted for whole mainlobe to rotate

% Calculate how many full mainlobe rotations will be detected in observationWindow
num_full_rotations_in_observation = floor(observationWindow / scanPeriod);
leftover_time_scans = mod(observationWindow, scanPeriod);

% Handle leftover time (time beyond full rotations)
if leftover_time_scans <= t_mainlobe_rotation
    mainlobe_leftover_time = leftover_time_scans;
else
    mainlobe_leftover_time = t_mainlobe_rotation;
end

% Calculate how many scans in total
if mainlobe_leftover_time <= t_mainlobe_rotation
    num_scans_in_observation = num_full_rotations_in_observation + 1;
else
    num_scans_in_observation = num_full_rotations_in_observation;
end

pulse_all_scans = pulse_detection_rotation*num_scans_in_observation;

% Calculate the time (indices as well) at which pulses occurs
% Index of scan_time array
%---------
% Time
pulses_mainlobe = zeros(2,pulse_all_scans);
pulses_mainlobe_amp = zeros(1,pulse_all_scans);
pulses_mainlobe_aoa = zeros(1,pulse_all_scans);

% loop through all scans
for scan_num = 1:num_scans_in_observation
    % Adjust first detection time based on scan number
    scan_shift = (scan_num-1)*pulse_detection_rotation;
    esm_first_detection = radar.TOA + (scan_num-1)*scanPeriod;
    time_shift = esm_first_detection - time_first_detection;
    % Shift scan_time array
    scan_time_shift = scan_time+time_shift;
    scan_time_mainlobe_shift = scan_time_mainlobe+time_shift;
    % Find all PRI indices and associated time over scan in mainlobe
    j = 2;
    pulses_mainlobe(1,1+scan_shift) = 1; % First pulse at first index
    pulses_mainlobe(2,1+scan_shift) = esm_first_detection; % First pulse at t = TOA
    for k = 2+scan_shift:pulse_detection_rotation + scan_shift
        % Calculate time for k-th pulse
        pulse_time = esm_first_detection +(j)*PRI;
        % Find closest index
        [~,closest_idx] = min(abs(scan_time_mainlobe_shift-pulse_time));
        pulses_mainlobe(1,k) = closest_idx;
        pulses_mainlobe(2,k) = scan_time_mainlobe_shift(closest_idx);
        j = j+1;
    end
    
    % Calculate associated Amp & AOA change based on TOA for each pulse
    m = 1;
    for i=1+scan_shift:pulse_detection_rotation + scan_shift
        pulse_time = esm_first_detection +(m)*PRI;
        [~,closest_idx] = min(abs(scan_time_mainlobe_shift-pulse_time));
        pulses_mainlobe_aoa(i) = mod(scanRate*pulses_mainlobe(2,i),360);
        pulses_mainlobe_amp(i) = amplitude_dB(closest_idx);
        m = m+1;
    end
end
% Generate PDWs for all pulses in single mainlobe
% PDW_mainlobe(pulse_detection_rotation) = struct();
for detected_pulse = 1:pulse_all_scans
    PDW_mainlobe(detected_pulse) = generatePDW(PDWParameters,pulses_mainlobe(2,detected_pulse),AOA+pulses_mainlobe_aoa(detected_pulse),pulses_mainlobe_amp(detected_pulse),radar.Freq,radar.PW);
end

%%


% Account for additional scans occuring in observation window
% Calculate how many full mainlobe rotations will be detected in observationWindow
num_full_rotations_in_observation = floor(observationWindow / scanPeriod);
leftover_time_scans = mod(observationWindow, scanPeriod);

% Handle leftover time (time beyond full rotations)
if leftover_time_scans <= t_mainlobe_rotation
    mainlobe_leftover_time = leftover_time_scans;
else
    mainlobe_leftover_time = t_mainlobe_rotation;
end

% Calculate how many scans in total
if leftover_time_scans <= t_mainlobe_rotation
    num_scans_in_observation = num_full_rotations_in_observation + 1;
else
    num_scans_in_observation = num_full_rotations_in_observation;
end

% Append detections for each additional scan
for scan_idx = 1:num_full_rotations_in_observation
    scan_offset_time = scan_idx*scanPeriod; % Time offset for each scan

    % Append detections for this scan
    for detected_pulse = 1:pulse_detection_rotation
        pulse_time_shifted = pulses_mainlobe(2,:) + scan_offset_time;
        % Find closest index after shifting
        [~, closest_idx] = min(abs(scan_time_mainlobe_shift - (pulses_mainlobe(2,1))));

        % Calculate AOA and amplitude for the shifted pulses
        shifted_aoa = mod(scanRate*pulse_time_shifted,360);
        shifted_amp(detected_pulse) = amplitude_dB(closest_idx);

        % Append new PDW for this scan's pulse
        PDW_mainlobe(end + 1) = generatePDW(PDWParameters, pulse_time_shifted(2,detected_pulse), ...
            AOA + shifted_aoa(detected_pulse), shifted_amp(detected_pulse), radar.Freq, radar.PW);
    end
end
%% Repeats due to scan
num_full_rotations_in_observation = floor(observationWindow/scanPeriod);
disp(['Number of full rotations: ',num2str(num_full_rotations_in_observation)]);
leftover_time_scans = mod(observationWindow,scanPeriod);
disp(['Leftover time: ',num2str(leftover_time_scans)]);
percentage_of_full_mainlobe = (t_mainlobe_rotation/scanPeriod)*100;
disp(['Mainlobe is: ', num2str(percentage_of_full_mainlobe),'% of full scan period']);
num_scans_in_observation = observationWindow/scanPeriod;
disp(['Number of rotations:',num2str(num_scans_in_observation)]);

if leftover_time_scans<=t_mainlobe_rotation
    mainlobe_leftover_time = leftover_time_scans;
else
    mainlobe_leftover_time = t_mainlobe_rotation;
end

after_mainlobe_leftover_time = observationWindow - (num_full_rotations_in_observation*scanPeriod + mainlobe_leftover_time);
disp(['Leftover time after next mainlobe: ',num2str(after_mainlobe_leftover_time)]);

if after_mainlobe_leftover_time >= 0
    num_scans_in_observation = num_full_rotations_in_observation + 1;
else
    num_scans_in_observation = num_full_rotations_in_observation; % Update to include any %
end

disp(num_scans_in_observation);