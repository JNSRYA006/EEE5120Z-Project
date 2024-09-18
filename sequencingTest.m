outputFileName = 'testSequencing.txt';

% Assuming delta_toa is your structure, extract pri_values
pri_values = cell(numel(delta_toa), 4);
for i = 1:numel(delta_toa)
    cluster = delta_toa(i);
    for j = 1:4
        pri = cluster.Levels(j).PRI;
        pri_values{i,j} = pri;
    end
end

% Call the sequencing function
sequenced_pulses = sequence_pulses(PDWImport, pri_values, radarParameters, outputFileName);



function sequenced_pulses = sequence_pulses(PDWImport, pri_values, radarParams, outputFileName)
    % Initialize the sequenced_pulses structure
    sequenced_pulses = struct([]);
    
    % Open the output file
    fid = fopen(outputFileName, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', outputFileName);
    end
    
    num_clusters = size(pri_values, 1);  % Number of clusters
    num_levels = size(pri_values, 2);    % Number of levels (should be 4)
    
    % Sort PDWImport by TOA for sequential processing
    [~, idx] = sort([PDWImport.TOA]);
    PDWImport = PDWImport(idx);
    
    % Process each cluster
    for cluster_idx = 1:num_clusters
        % Extract PRI candidates for the cluster
        pri_candidates = zeros(num_levels, 3);  % Initialize the PRI matrix
        for level = 1:num_levels
            pri_value = cell2mat(pri_values(cluster_idx, level));
            pri_candidates(level, :) = pri_value / level;  % Store PRI candidates for this level
        end
        
        % Compute the intersection of PRI candidates across levels
        pri_candidate = pri_candidates(1, :);
        for level = 2:num_levels
            pri_candidate = intersect(pri_candidate, pri_candidates(level, :));
        end
        
        % Start the sequencing process for the current cluster
        pdw_idx = 1;  % Start from the first PDW
        scan_start_TOA = PDWImport(1).TOA;  % Initial TOA
        
        while pdw_idx <= length(PDWImport) && PDWImport(pdw_idx).TOA < scan_start_TOA + radarParams.radar1.scanPeriod
            % Start with the first pulse
            first_TOA = PDWImport(pdw_idx).TOA;
            current_sequence = PDWImport(pdw_idx);  % Initialize the current sequence
            
            pdw_idx = pdw_idx + 1;  % Move to the next PDW
            
            % Look for pulses based on PRI candidates
            max_missing_pulses = 5;
            tolerance = 5e-7;
            pulse_found = true;
            sequence_length = 1;
            missing_pulses_count = 0;
            
            while pulse_found && missing_pulses_count <= max_missing_pulses
                expected_TOA = first_TOA + pri_candidate * (sequence_length + missing_pulses_count);
                
                % Get the TOA values of all PDWs
                pdw_toa_values = [PDWImport.TOA];
                
                % Ensure element-wise subtraction
                next_pulse_idx = find(abs(pdw_toa_values - expected_TOA) < tolerance);
                
                if ~isempty(next_pulse_idx)
                    % Add the found pulse to the sequence
                    current_sequence(end+1) = PDWImport(next_pulse_idx); %#ok<AGROW>
                    
                    % Mark the pulse as used by setting TOA to Inf
                    PDWImport(next_pulse_idx).TOA = Inf;
                    
                    % Reset missing pulses count and increase sequence length
                    missing_pulses_count = 0;
                    sequence_length = sequence_length + 1;
                else
                    % No pulse found, increment missing pulses count
                    missing_pulses_count = missing_pulses_count + 1;
                end
            end
            
            % If valid sequence found, add it to sequenced_pulses
            if ~isempty(current_sequence)
                sequenced_pulses = [sequenced_pulses, current_sequence]; %#ok<AGROW>
            end
            
            % Break if the expected TOA exceeds the scan window
            if expected_TOA > scan_start_TOA + radarParams.radar1.scanPeriod
                break;
            end
        end
    end
    
    % Close the output file
    fclose(fid);
end
