function sequencing(PDWs,delta_toa,radarParams,outputFileName)
   
    % Extract PRI values from PRICandidates structure
    pri_values = cell(numel(delta_toa),4);
    for i = 1:numel(delta_toa)
        cluster = delta_toa(i);
        for j = 1:4
            pri = cluster.Levels(j).PRI;
            pri_values{i,j} = pri;
        end
    end

    % Open the output file
    fid = fopen(outputFileName, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', outputFileName);
    end

    num_radars = length(radarParams);
    % Go through each cluster and identify radars based on PRI candidates
    % determined in the differential TOA process
    % Define the number of clusters and levels
    num_cluster = length(clusters);
    num_levels = 4;
    
    % Initialize the matrix for PRI candidates
    % Assuming each PRI candidate is a row vector of length 3
    pri_candidates = zeros(num_levels, 3);
    % pri_selected_candidates = zeros(4,1);
    
    % Sort PDWImport by TOA for sequential processing
    [~, idx] = sort([PDWs.TOA]);
    PDWs = PDWs(idx); % Sort PDWImport by TOA

    % pulses_to_remove = [];

    % Loop through each cluster
    for cluster_idx = 1:num_cluster
        % Loop through each level
        for level = 1:num_levels
            % Extract PRI candidates at each level
            pri_value = cell2mat(pri_values(cluster_idx, level));
            % Compute the pri_candidate
            pri_candidate = pri_value/level;
            % Determine the row index in the matrix
            row_idx = (cluster_idx - 1) * num_levels + level;
            pri_candidates(level, :) = pri_candidate;
        end
        % Account for precision issues
        precision = 9;
        pri_candidates = round(pri_candidates,precision);
        % Determine the PRI value which is common throughout all 4 levels
        pri_candidate = pri_candidates(1,:);
        for row_idx = 2:size(pri_candidates,1)
            pri_candidate = intersect(pri_candidate,pri_candidates(row_idx,:));
        end
        pri_selected_candidates = [];
        pri_selected_candidates = [pri_selected_candidates,pri_candidate];
        pri_candidates = zeros(num_levels, 3);

        % Go through merged list (PDWs) and start at first TOA. Look at a PRI
        % candidate away and see if a pulse exists. If it does, add to
        % seuqenced list, then go forward until no pulse exists (do 5
        % more to mitigate pulse loss and uncertainty).
        % After succesfully finding PRIs for single scan, discard these
        % from the PDWImport and begin processing again from next first
        % TOA value.
        
        % Initialize sequenced pulses list
        sequenced_pulses = struct([]);
        
        % Initialize a flag for when to stop within the scan
        scan_start_TOA = PDWs(1).TOA;

        % Keep track of pulses to remove
        pulses_to_remove = [];
        
        % Process only pulses within the single scan (scanPeriod)
        pdw_idx = 1; % Index for PDWs
        while pdw_idx <= length(PDWs) && PDWs(pdw_idx).TOA < scan_start_TOA + radarParams.radar3.scanPeriod
            % Start from the first pulse (TOA) in PDWs
            first_TOA = PDWs(pdw_idx).TOA;
            
            % Initialize a temporary list to hold the current sequence of pulses
            current_sequence = PDWs(pdw_idx); % Start with the first pulse
            
            % Move to the next PDW in PDWs
            pdw_idx = pdw_idx + 1;
            
            % Loop over each PRI candidate
            for pri_idx = 1:numel(pri_selected_candidates)
                % Get the current PRI candidate to test
                PRI = pri_selected_candidates(pri_idx);
                
                % Start looking for pulses based on the PRI candidate
                missing_pulses_count = 0;
                pulse_found = true;
                sequence_length = 1; % Start with one pulse in the sequence
                max_missing_pulses = 5; % Set the maximum number of missing pulses allowed
                tolerance = PDWParameters.TOA.resolution;
                while pulse_found && missing_pulses_count <= max_missing_pulses
                    % Calculate expected TOA based on PRI and current sequence length
                    expected_TOA = first_TOA + PRI * (sequence_length+missing_pulses_count);
                    
                    % Find the pulse that matches the expected TOA within a
                    % tolerance of uncertainty
                    TOA_diffs = abs([PDWs.TOA] - expected_TOA);
                    next_pulse_idx = find(TOA_diffs < tolerance, 1);
                    % next_pulse_idx = find(abs([PDWs.TOA] - expected_TOA) < tolerance,1);
                    
                    if ~isempty(next_pulse_idx)
                        % Add the found pulse to the sequence
                        current_sequence(end+1) = PDWs(next_pulse_idx);
                        
                        % Mark the found pulse to remove
                        pulses_to_remove(end+1) = next_pulse_idx;
                        
                        % Reset missing pulses count and increment sequence length
                        missing_pulses_count = 0;
                        sequence_length = sequence_length + 1;
                    else
                        % If no pulse found, increment missing pulses count
                        missing_pulses_count = missing_pulses_count + 1;
                        % Exit if no pulses found
                        if missing_pulses_count > max_missing_pulses
                            pulse_found = false;
                        end
                        
                    end
                end
            end
            
            % If a valid sequence is found, add it to the sequenced pulses list
            if ~isempty(current_sequence)
                sequenced_pulses = [sequenced_pulses, current_sequence]; %#ok<AGROW>
            end
            
            % Break the loop if the TOA exceeds the current scan window
            if expected_TOA > scan_start_TOA + radarParams.radar3.scanPeriod
                break;
            end
        end
        % Remove pulses
        PDWs(pulses_to_remove) = [];

        % Reset
        pulses_to_remove = [];
    end
    


end