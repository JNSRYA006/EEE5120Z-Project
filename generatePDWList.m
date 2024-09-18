function PDW = generatePDWList(PDWParameters,radarParams,radarName,observationWindow,pulseLossPercentage)

    lockOnTx = 200e-3; % ms
    sectorScanWidth = 10; % deg

    % Check that scan mode is of correct type
    validScanModes = ["none","circular","sector","lock"];
    % Check that correct radar is requested for generating PDWs
    validRadars = ["radar1","radar2","radar3","radar4"];

    % Extract scanMode from all radari fields in radarParameters
    % radarCells = struct2cell(radarParams);
    radarCells = struct2cell(radarParams);
    % radars = string(cellfun(@(x) x.Field, radarParams, 'UniformOutput', false));
    scanModes = string(cellfun(@(x) x.scanMode, radarCells, 'UniformOutput', false));

    if ~ismember(lower(scanModes),validScanModes)
        error('Invalid scanMode. Allowed values are: "none", "circular", "sector", "lock" [Non case-sensitive]');
    end

    % Extract radar
    if ismember(radarName, validRadars)
        radar = radarParams.(radarName);
    else
        error('Field "%s" does not exist in the struct.', radarName);
    end

    % Check that scan mode is of correct type
    validScanModes = ["none","circular","sector","lock"];

    % Extract scanMode from all radari fields in radarParameters
    radarCells = struct2cell(radarParams);
    scanModes = string(cellfun(@(x) x.scanMode, radarCells, 'UniformOutput', false));

    if ~ismember(lower(scanModes),validScanModes)
        error('Invalid scanMode. Allowed values are: "none", "circular", "sector", "lock" [Non case-sensitive]');
    end

    switch lower(radar.scanMode)
        case "circular"
            PDW = circularScan(PDWParameters,radar,observationWindow);
        case "sector"
            PDW = sectorScan(PDWParameters,radar,observationWindow,sectorScanWidth);
        case "lock"
            PDW = lockOnScan(PDWParameters,radar,lockOnTx);
        otherwise
            error('Invalid scanMode. Allowed values are: "none", "circular", "sector", "lock" [Non case-sensitive]');
    end
    % 2.6 - Add pulse loss to PDW lists
    PDW = applyPulseLoss(PDW,pulseLossPercentage);

end

function PDW = circularScan(PDWParameters,radar,observationWindow)
    % Extract radar parameters
    PRI = radar.PRI;
    beamwidth = radar.Beamwidth; 
    mainlobe = 2*beamwidth;
    scanPeriod = radar.scanPeriod; % s
    scanRate = 360/scanPeriod;
    sensitivity = -50; % dBm
    
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
    scan_time_shift = scan_time + time_shift;
    scan_time_mainlobe_shift = scan_time_mainlobe + time_shift;
    
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
        
        % Calculate associated Amp change based on TOA for each pulse
        m = 1;
        for i=1+scan_shift:pulse_detection_rotation + scan_shift
            pulse_time = esm_first_detection +(m)*PRI;
            [~,closest_idx] = min(abs(scan_time_mainlobe_shift-pulse_time));
            pulses_mainlobe_amp(i) = amplitude_dB(closest_idx);
            m = m+1;
        end
    end
    % Generate PDWs for all pulses in single mainlobe
    % PDW_mainlobe(pulse_detection_rotation) = struct();
    for detected_pulse = 1:pulse_all_scans
        PDW_mainlobe(detected_pulse) = generatePDW(PDWParameters,pulses_mainlobe(2,detected_pulse),radar.AOA,pulses_mainlobe_amp(detected_pulse),radar.Freq,radar.PW);
    end
    PDW = PDW_mainlobe;

end

function PDW = sectorScan(PDWParameters,radar,observationWindow,scanWidth)

    % Extract radar parameters
    first_intercept = radar.TOA;
    PRI = radar.PRI;
    beamwidth = radar.Beamwidth;
    mainlobe = 2*beamwidth;
    scanPeriod = radar.scanPeriod; % s
    scanRate = scanWidth/scanPeriod;
    sensitivity = -50; % dBm
    
    % Num of pulses per beamwidth rotation
    disp('-----------------------------------------------');
    fprintf('Scan rate of Tx radar: %.1f deg/s\n',scanRate);
    t_beamwidth_rotation = beamwidth/scanRate; % Time for 3dB beamwidth to rotate through
    t_mainlobe_rotation = mainlobe/scanRate; % Time for mainlobe to rotate through full size
    fprintf('%.3f ms for rotation through %.2f deg (mainlobe)\n',t_mainlobe_rotation*1e3,mainlobe);
    pulse_beamwidth_rotation = floor(t_beamwidth_rotation/PRI); % Number of pulses emitted for 3dB beamwidth to rotate
    pulse_mainlobe_rotation = floor(t_mainlobe_rotation/PRI); % Number of pulses emitted for whole mainlobe to rotate
    fprintf('%.0f full pulses emitted in rotation through full mainlobe\n',pulse_mainlobe_rotation);
    disp('-----------------------------------------------');
    
    pulses = floor((observationWindow-first_intercept)/PRI);
    pdw_idx = 1;
    
    for pulse = 1:pulses
    
        TOA = (pulse-1)*PRI + first_intercept;
        theta_shift = scanRate*TOA;
        theta_sector = mod(theta_shift+(9 - beamwidth),2*scanWidth);
    
        % Compute the function f(theta) based on AOA
        amplitude = sin(((theta_sector - beamwidth+9)*pi)/(2*beamwidth));
        amplitude_dB = 20*log10(amplitude);
        amplitude_dB = radar.peakAmp + amplitude_dB;
    
        % Check if AOA is within the detection sector and above sensitivity threshold
        esm_detection_angle_1 = 9;
        esm_detection_angle_2 = 11;
    
        if (abs(theta_sector - esm_detection_angle_1) < beamwidth) || (abs(theta_sector - esm_detection_angle_2) < beamwidth)
            if amplitude_dB >= sensitivity
                pulse_amplitude = amplitude_dB;
            else
                pulse_amplitude = -Inf;
            end
        else
            pulse_amplitude = -Inf;
        end
        % Only generate PDW for values > sensitivity
        if pulse_amplitude ~= -Inf
            % Generate PDW
            PDW_mainlobe(pdw_idx) = generatePDW(PDWParameters,TOA,radar.AOA,pulse_amplitude,radar.Freq,radar.PW);
            pdw_idx = pdw_idx + 1;
        end
    end
    PDW = PDW_mainlobe;
end


function PDW = lockOnScan(PDWParameters,radar,transmissionTime)

    %% Initialise parameters    
    AOA = radar.AOA;
    PRI = radar.PRI;
    peakAmp = radar.peakAmp;
    
    %% Calculate Tx Mainlobe
    % Since this is lockon, the peak amplitude will be detected every PRI
    % after the first intercept time at the designated AOA.
    % The only parameter that will change is the TOA.
    
    % Need to offset this first detection to occur at t = TOA (Based on Table
    % 2)
    % Fist detect per Table 2
    esm_first_detection = radar.TOA;
    
    t_pulses_detected = transmissionTime - esm_first_detection;
    % disp(t_pulses_detected);
    
    % Calculate the number of pulses detected in transmission time
    t = esm_first_detection;
    pulse_detection_transmitting = 1; % First pulse detected at t = first intercept
    pulses_mainlobe_toa = zeros(1,floor(t_pulses_detected/PRI));
    while t <= transmissionTime
        pulses_mainlobe_toa(pulse_detection_transmitting) = t;
        pulse_detection_transmitting = pulse_detection_transmitting + 1;
        t = t + PRI;
    end
    
    % Account for last pulse being added but not detected
    pulse_detection_transmitting = pulse_detection_transmitting - 1;
    
    disp('-----------------------------------------------');
    fprintf('%.0f full pulses emitted during transmission time\n',pulse_detection_transmitting);
    disp('-----------------------------------------------');
    
    %% Generate PDWs
    pulses_mainlobe_amp = ones(1,pulse_detection_transmitting)*peakAmp;
    pulses_mainlobe_aoa = ones(1,pulse_detection_transmitting)*AOA;
    
    for detected_pulse = 1:pulse_detection_transmitting
        PDW_mainlobe(detected_pulse) = generatePDW(PDWParameters,pulses_mainlobe_toa(detected_pulse),pulses_mainlobe_aoa(detected_pulse),pulses_mainlobe_amp(detected_pulse),radar.Freq,radar.PW);
    end
    PDW = PDW_mainlobe;

end


function PDW_w_losses = applyPulseLoss(PDW, lossPercentage)

    % Check lossPercentage is within allowable bounds
    if lossPercentage < 0 || lossPercentage > 100
        error('Loss percentage must be between 0 and 100.');
    end

    % Calculate the number of pulses to get rid of
    num_pulses = length(PDW);
    num_pulses_remove = round(num_pulses*(lossPercentage/100));
    % Select the indices to remove
    indices_remove = randperm(num_pulses,num_pulses_remove);

    % Remove the selected pulses
    PDW_w_losses = PDW;
    PDW_w_losses(indices_remove) = [];

    % Output verification
    fprintf('Removed %d of %d pulses.\n',num_pulses_remove,num_pulses);
end