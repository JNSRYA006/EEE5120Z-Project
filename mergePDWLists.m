function mergedPDW = mergePDWLists2(varargin)
    % mergePDWLists - Merges PDW lists, sorts them by TOA, and handles overlapping pulses.

    % Concatenate all PDW structs from different lists into one array
    merged_pdw = [varargin{:}];  % More efficient than concatenating inside the loop

    % Sort the PDWs by TOA
    [~, idx] = sort([merged_pdw.TOA]);
    merged_pdw = merged_pdw(idx);

    % Preallocate for result_pdw (assuming worst case, no overlaps)
    result_pdw = repmat(struct('TOA', [], 'AOA', [], 'Amp', [], 'Freq', [], 'PW', []), size(merged_pdw));

    % Start merging the PDWs
    current_pdw = merged_pdw(1);
    result_index = 1;  % Index to keep track of where to insert in result_pdw

    for i = 2:length(merged_pdw)
        next_pdw = merged_pdw(i);
        
        % Check if the current pulse overlaps with the next one
        if next_pdw.TOA < current_pdw.TOA + current_pdw.PW
            % Merge pulses (overlapping case)
            current_pdw.PW = max(current_pdw.TOA + current_pdw.PW, next_pdw.TOA + next_pdw.PW) - current_pdw.TOA;
            
            % Zero out other fields (AOA, Amp, Freq)
            current_pdw.AOA = 0;
            current_pdw.Amp = 0;
            current_pdw.Freq = 0;
        else
            % No overlap: save the current PDW and move to the next one
            result_pdw(result_index) = current_pdw;
            result_index = result_index + 1;
            current_pdw = next_pdw;
        end
    end
    
    % Append the final PDW
    result_pdw(result_index) = current_pdw;
    
    % Trim the preallocated result_pdw to remove unused entries
    mergedPDW = result_pdw(1:result_index);
end
