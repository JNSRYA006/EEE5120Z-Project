function differences = differentialTOAAnalysis(clusters,barWidth,plotHisto)
    % Differential TOA analysis with CDIF histograms for each cluster
    % Inputs:
    %   clusters: Cell array of clustered PDWs (result from the clustering step)
    %   barWidth: Width of histogram bars for ΔTOA histograms
    %   plotHisto: Whether or not to plot of CDIF histograms
    %
    % This function constructs CDIF histograms up to the 4th level for each cluster.
    
    % Loop through each cluster to perform ΔTOA analysis
    differences(length(clusters)) = struct;
    for c = 1:length(clusters)
        pdw_cluster = clusters{c};
        num_pdw = length(pdw_cluster);
        
        if num_pdw < 4
            fprintf('Cluster %d has fewer than 4 PDWs, skipping...\n', c);
            continue;
        end
        
        % Extract TOA from the PDWs in this cluster
        TOAs = [pdw_cluster.TOA];
        % Define the number of levels to investigate
        levels = 4;

        % Initialise a structure for the current cluster
        differences(c).Cluster = c;
        % differences(c).level = struct();

        % Loop through each spacing level
        for level = 1:levels
            % Initialize an array to store the differences for the current spacing
            current_delta = zeros(1, length(TOAs) - level);
            
            % Compute differences for the current spacing
            for i = 1:(length(TOAs) - level)
                current_delta(i) = TOAs(i + level) - TOAs(i);
            end
            
            % Store the differences in the structure array
            % Get most common deltaTOA value and store (with uncertainty) as
            % a PRI candidate. Returned for sequencing
            differences(c).Levels(level).delta = current_delta;
            bin_edges = min(current_delta):barWidth:max(current_delta);
            h = histogram(current_delta, 'BinWidth', barWidth,'BinEdges',bin_edges);
            % h = histogram(current_delta, 'BinWidth', barWidth);

            % Get bin counts
            counts = h.BinCounts;
            % Find the index of the bin with highest frequency
            max_idx = find(counts == max(counts));
            % Use bin edges to determine bounds of highest f bin
            PRI_vals = [bin_edges(max_idx) bin_edges(max_idx+1)];
            % Take best guess as the middle value of the bin
            best_guess = (PRI_vals(1)+PRI_vals(2))/2;
            % Write best guess PRI to array with uncertainty above and
            % below based on resolution of bars
            differences(c).Levels(level).PRI = [best_guess-barWidth/2 best_guess best_guess+barWidth/2];
            
            fprintf('Completed ΔTOA analysis for Cluster %d, Level %d.\n',c,level);
            
            if plotHisto == 1
                % Create array to define bin edges
                % bin_edges = min(current_delta):barWidth:max(current_delta);
                % fprintf('Creating CDIF histogram for Cluster %d...\n', c);
                figureHandle = figure('Name', sprintf('Histogram for Cluster %d Level %d', c, level),'NumberTitle', 'off', 'Visible', 'on');
                clf(figureHandle);
                histogram(current_delta, 'BinWidth', barWidth,'BinEdges',bin_edges);
                title(sprintf('CDIF Histogram (Level %d) for Cluster %d', level,c));
                xlabel('ΔTOA [s]');
                ylabel('Number of Pulses');
                grid on;
                ax = gca;
                ax.XAxis.Exponent = -6;
                axis auto;
                drawnow;
                figureName = sprintf('/home/ryan/Documents/EEE5120Z/Project/Figures/Deint/cluster_%d_cdif_level_%d.tex',c,level);
                matlab2tikz(figureName);
                fprintf('Figure: %s, saved in Figures/ folder\n',figureName);
                pause(1);
            end
        end
    end
end