function PRICandidates = selectPRICandidates(deltaTOA)
    % Inputs:
    % delta_toa - 1xn structure containing the cluster information and delta levels
    % Output:
    % PRI_candidates - Cell array where each cell contains the best PRI candidate for each cluster

    num_clusters = length(deltaTOA);
    PRICandidates = cell(1, num_clusters); % Preallocate cell array for PRI candidates

 % Loop through each cluster
    for cluster_idx = 2:num_clusters
        cluster_data = deltaTOA(cluster_idx); % Get data for the current cluster
        
        % Collect all delta values across the 1x4 Levels structure
        delta_values = [];
        for level_idx = 1:length(cluster_data.Levels)
            delta_values = [delta_values, cluster_data.Levels(level_idx).delta(level_idx)];
        end
        
        % Find the most frequent delta value
        [unique_deltas, ~, idx] = unique(delta_values); % Get unique delta values and their indices
        freq_counts = histcounts(idx, 1:numel(unique_deltas)+1); % Count occurrences of each unique delta
        
        % Find the delta value that occurs the most
        [~, max_idx] = max(freq_counts); % Get the index of the most frequent delta value
        most_frequent_delta = unique_deltas(max_idx); % Retrieve the most frequent delta value
        
        % Select the level that contains the most frequent delta value
        best_level_idx = find([cluster_data.Levels.delta] == most_frequent_delta, 1);
        
        % Store the best PRI candidate (based on the most frequent delta) for this cluster
        PRI_candidates{cluster_idx} = best_level_idx;
    end
end