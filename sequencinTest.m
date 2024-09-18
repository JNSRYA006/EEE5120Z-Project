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



function sequenced_pulses = sequence_pulses(PDWs, pri_values, radarParams, outputFileName)
    % Extract PRI values from PRICandidates structure
    pri_candidates = zeros(numel(pri_values), 4, 3);  % Assuming 4 levels and 3 PRI values per level
    for i = 1:numel(pri_values)
        for j = 1:4
            pri_candidates(i, j, :) = pri_values{i, j} / j;  % Division by level may need to be adjusted
        end
    end

    % Initialize sequenced pulses list
    sequenced_pulses = struct([]);
    
    % Sort PDWImport by TOA for sequential processing
    [~, idx] = sort([PDWs.TOA]);
    PDWs = PDWs(idx); % Sort PDWImport by TOA
    
    % Set tolerance for TOA matching
    tolerance = radarParams.TOA.resolution * 5;  % Increase tolerance for more flexibility
    
    % Process pulses within the single scan (scanPeriod)
    pdw_idx = 1;
    while pdw_idx <= length(PDWs)
        first_TOA = PDWs(pdw_idx).TOA;
        current_sequence = PDWs(pdw_idx); % Initialize sequence with first pulse
        pdw_idx = pdw_idx + 1;
        
        % Try each PRI candidate to find matching pulses
        for pri_level = 1:4
            PRI = pri_candidates(:, pri_level, :);
            sequence_length = 1;
            missing_pulses = 0;
            max_missing_pulses = 5;
            pulse_found = true;
            
            while pulse_found && missing_pulses <= max_missing_pulses
                expected_TOA = first_TOA + PRI * sequence_length;  % Expected next pulse's TOA
                
                % Find next pulse within the tolerance window
                next_pulse_idx = find(abs([PDWs.TOA] - expected_TOA) < tolerance, 1);
                
                if ~isempty(next_pulse_idx)
                    current_sequence(end+1) = PDWs(next_pulse_idx); %#ok<AGROW>
                    % Mark pulse as processed
                    PDWs(next_pulse_idx).processed = true;
                    sequence_length = sequence_length + 1;
                    missing_pulses = 0;  % Reset missing pulse count
                else
                    missing_pulses = missing_pulses + 1;
                end
            end
        end
        
        % If a valid sequence is found, add it to the sequenced pulses list
        if ~isempty(current_sequence)
            sequenced_pulses = [sequenced_pulses, current_sequence]; %#ok<AGROW>
        end
    end

    % Remove processed pulses from PDWImport
    PDWs = PDWs(~[PDWs.processed]);
    
    % Save sequenced pulses to file
    save_sequenced_pulses(outputFileName, sequenced_pulses);
end

function save_sequenced_pulses(outputFileName, sequenced_pulses)
    fid = fopen(outputFileName, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', outputFileName);
    end
    
    for i = 1:length(sequenced_pulses)
        fprintf(fid, 'Pulse: %d, TOA: %f\n', i, sequenced_pulses(i).TOA);
    end
    fclose(fid);
end
