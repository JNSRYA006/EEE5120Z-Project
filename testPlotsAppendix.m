radar = radarParameters.radar2;

% Define the beamwidth in degrees (3-dB beamwidth)
beamwidth = radar.Beamwidth; 
scanPeriod = radar.scanPeriod; % s
scanRate = 360/scanPeriod;
fprintf('%.0f deg/s\n',scanRate);
sensitivity = -50; % dB = -50 dBm

% Num of pulses per beamwidth rotation
t_beamwidth_rotation = beamwidth/scanRate;
fprintf('%.3f ms for rotation through %.2f deg (beamwdith)\n',t_beamwidth_rotation*1e3,beamwidth);
PRI = radar.PRI;
pulse_beamwidth_rotation = t_beamwidth_rotation/PRI;
fprintf('%.0f full pulses emitted in rotation through beamwidth\n',ceil(pulse_beamwidth_rotation));
disp('-----------------------------------------------')


n = 1000000;
% Define the angle theta from 0 to 360 degrees
theta = linspace(0, 2*beamwidth, n);
theta_rad = deg2rad(theta);
theta_mainlobe = theta+radar.AOA;
theta_mainlobe(theta_mainlobe>2*beamwidth+radar.AOA) = 0;

% Compute the function f(theta)
% detection_angle_shift = radar.TOA*scanRate;
% time_shift = detection_angle_shift - theta(detection_index_first);
amplitude = sin((theta*pi)/(2*beamwidth));

% Compute the amplitude in dB
amplitude_dB = 20*log10(amplitude)+(radar.peakAmp);

% Find first index at which amplitude crosses sensitivity
% detection_index_first = find(amplitude_dB > sensitivity & theta >= (radar.AOA-1) & theta <= (radar.AOA+2*beamwidth+1),1,"first");
detection_index_first = find(amplitude_dB > sensitivity,1,"first");
detection_index_first_offset = find(amplitude_dB > sensitivity,1,"first");
detection_index_last = find(amplitude_dB > sensitivity,1,"last");
% detection_index_last = find(amplitude_dB > sensitivity & theta >= (radar.AOA-1) & theta <= (radar.AOA+2*beamwidth+1),1,"last");

fprintf('First detection > sensitivity occurs at index: %.0f. Amp [dB] = %.4f, Thetha [deg] = %.4f.\n',detection_index_first,amplitude_dB(detection_index_first),theta(detection_index_first));
% disp(detection_index_first);
% disp(amplitude_dB(detection_index_first));
% disp(theta(detection_index_first));

fprintf('Last detection > sensitivity occurs at index: %.0f. Amp [dB] = %.4f, Thetha [deg] = %.4f.\n',detection_index_last,amplitude_dB(detection_index_last),theta(detection_index_last));
% disp(detection_index_last);
% disp(amplitude_dB(detection_index_last));

% Shifting to first intercept
detection_angle_shift = radar.TOA*scanRate;
time_shift = -abs(detection_angle_shift + theta(detection_index_first_offset));
thetaShifted = theta + time_shift;
thetaMainLobeShifted = theta_mainlobe + time_shift;
fprintf('Angle shift due to TOA from scan rate = %.4f.\n',detection_angle_shift);
fprintf('Angle shift due to TOA from theta index = %.4f.\n',theta(detection_index_first));


amplitude_dB_shifted = amplitude_dB + time_shift;

aoa_lims = [radar.AOA 2*beamwidth+radar.AOA+2*time_shift];
aoa_range = aoa_lims(2)-aoa_lims(1);
fprintf('Range of AOA detected: %.4f deg\n',aoa_range);

% Plot the amplitude in dB
figure;
plot(thetaShifted, amplitude_dB, 'LineWidth', 2);
xlabel('Angle \theta (degrees)');
ylabel('Amplitude (dB)');
title(['Amplitude of Radar Mainlobe with ' num2str(beamwidth) '° Beamwidth']);
grid on;
% xlim([radar.AOA 2*beamwidth+radar.AOA+2*time_shift]);
% xlim([0 2*beamwidth-detection_angle_shift]);
% ylim([-100 max(amplitude_dB)]); % Adjust this based on the expected range of amplitudes
hold on;
yline(sensitivity,'r--','LineWidth',1,'Label',['Sensitivity:' num2str(sensitivity) 'dB']);
xline(beamwidth,'LineStyle','--','LineWidth',1,'Label','Boresight');

%% Convert to scan period
boresight_time = (beamwidth/360)*scanPeriod;
scanTime = linspace(0,scanPeriod,n);
scanTimeBeamwidth = linspace(0,2*boresight_time,n);
PRI = radar.PRI;
t = radar.TOA;
disp(scanTimeBeamwidth(detection_index_last));

% Shifting to first intercept
scan_time_shift = radar.TOA - scanTime(detection_index_first);
scanTimeShifted = scanTime + scan_time_shift;
boresight_time_shift = boresight_time + scan_time_shift;

% Calculating time differences of interest
t_first = scanTimeBeamwidth(detection_index_first);
disp(t_first);
t_detect = 2*boresight_time - t_first;
disp(t_detect);
first_intercept = radar.TOA; % Actually t_first based on Table 2, because that is first sensitivity crossing
last_intercept = first_intercept + t_detect;
fprintf('First intercept: %.3f; Last intercept: %.3f (ms)\n',first_intercept*1e3,last_intercept*1e3);
%% 

% Plot the amplitude in dB
figure;
% plot(scanTime*1e3, amplitude_dB, 'LineWidth', 2);
% hold on;
plot(scanTimeShifted*1e3, amplitude_dB,'LineWidth', 1);
% hold off;
xlabel('Time (ms)');
ylabel('Amplitude (dB)');
title(['Amplitude of Radar Mainlobe with ' num2str(beamwidth) '° Beamwidth']);
grid on;
xlim([0 (2*boresight_time)*1e3]);
% ylim([-100 max(amplitude_dB)]); % Adjust this based on the expected range of amplitudes
hold on;
yline(sensitivity,'r--','LineWidth',1,'Label',['Sensitivity:' num2str(sensitivity) 'dB']);
xline(boresight_time*1e3,'LineStyle','--','LineWidth',1,'Label','Boresight');
% xline(radar.TOA*1e3,'g--','LineWidth',1,'Label','First intercept');
xline(scanTimeBeamwidth(detection_index_first)*1e3,'g--','LineWidth',1,'Label','First intercept');
disp(t*1e3);
disp(PRI*1e3);
while t <= t_detect
    xline(t*1e3,'b--','LineWidth',0.01,'Alpha',0.5);
    disp(t);
    t = t+PRI;
end
