%% Initialise parameters
radar = radarParameters.radar2;

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
%% Generate PDWs
% Calculate the time (indices as well) at which pulses occurs
% Index of scan_time array
%---------
% Time
pulses_mainlobe = zeros(2,pulse_detection_rotation);
pulses_mainlobe_amp = zeros(1,pulse_detection_rotation);
pulses_mainlobe_aoa = zeros(1,pulse_detection_rotation);
pulses_mainlobe(1,1) = 1; % First pulse at first index
pulses_mainlobe(2,1) = esm_first_detection; % First pulse at t = TOA
% Find all PRI indices and associated time over scan in mainlobe
for k = 2:pulse_detection_rotation
    % Calculate time for k-th pulse
    pulse_time = esm_first_detection +(k-1)*PRI;
    % Find closest index
    [~,closest_idx] = min(abs(scan_time_mainlobe_shift-pulse_time));
    pulses_mainlobe(1,k) = closest_idx;
    pulses_mainlobe(2,k) = scan_time_mainlobe_shift(closest_idx);
end

% Calculate associated Amp & AOA change based on TOA for each pulse
for i=1:pulse_detection_rotation
    pulse_time = esm_first_detection +(i)*PRI;
    [~,closest_idx] = min(abs(scan_time_mainlobe_shift-pulse_time));
    pulses_mainlobe_aoa(i) = scanRate*pulses_mainlobe(2,i);
    pulses_mainlobe_amp(i) = amplitude_dB(closest_idx);
end

% Generate PDWs for all pulses in single mainlobe
% PDW_mainlobe(pulse_detection_rotation) = struct();
for detected_pulse = 1:pulse_detection_rotation
    PDW_mainlobe(detected_pulse) = generatePDW(PDWParameters,pulses_mainlobe(2,detected_pulse),AOA+pulses_mainlobe_aoa(detected_pulse),pulses_mainlobe_amp(detected_pulse),radar.Freq,radar.PW);
end

%% Plot over scan time
% Plot the amplitude in dB
figure;
plot(scan_time_mainlobe_shift*1e3, amplitude_dB,'b--','LineWidth', 1);
hold on;
plot(scan_time_mainlobe*1e3, amplitude_dB,'k','LineWidth', 1);
xlabel('Time [ms]');
ylabel('Amplitude [dB]');
title(['Amplitude of Radar Mainlobe with ' num2str(beamwidth) '° Beamwidth']);
grid on;
% xlim([5*time_shift t_mainlobe_rotation*1e3]);
ylim([-80 max(amplitude_dB)]); % Adjust this based on the expected range of amplitudes
hold on;
yline(sensitivity,'r--','LineWidth',1,'Label',['Sensitivity:' num2str(sensitivity) 'dB']);
xline(t_boresight*1e3,'LineStyle','--','LineWidth',1,'Label','Boresight');
xline(esm_first_detection*1e3,'g--','LineWidth',1,'Label','First intercept (Table 2)');
% xline(scanTimeBeamwidth(detection_index_first)*1e3,'g--','LineWidth',1,'Label','First intercept');

%% Plot over angle
figure;
plot(theta_mainlobe_above_sensitivity, amplitude_dB_above_sensitivity, 'LineWidth', 1);
xlabel('Beamwidth, \theta [deg]','Interpreter','tex');
ylabel('Amplitude [dB]');
title(['Amplitude of Radar Mainlobe with ' num2str(beamwidth) '° Beamwidth']);
grid on;
% xlim([radar.AOA 2*beamwidth+radar.AOA+2*time_shift]);
% xlim([0 2*beamwidth-detection_angle_shift]);
ylim([-80 max(amplitude_dB)]); % Adjust this based on the expected range of amplitudes
hold on;
yline(sensitivity,'r--','LineWidth',1,'Label',['Sensitivity:' num2str(sensitivity) 'dB']);
xline(beamwidth,'LineStyle','--','LineWidth',1,'Label','Boresight','LabelOrientation','horizontal','LabelVerticalAlignment','bottom');