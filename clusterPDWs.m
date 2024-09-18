function clusters = clusterPDWs(PDWParams,PDWStruct,scatterPlot)
    % Parameters for clustering
    time_slice = 100e-3;  % 100 ms
    
    % Extract the required fields: AOA and Freq for clustering
    AOA = [PDWStruct.AOA];
    Freq = [PDWStruct.Freq];
    TOA = [PDWStruct.TOA];  % TOA is used for time slicing, not clustering
    
    % Apply time slice: filter PDWs within the first 100 ms
    valid_pdw = TOA <= time_slice;
    PDWStruct = PDWStruct(valid_pdw);
    
    % Cluster initialization
    clusters = {};
    cell_size = [];

    % Initial cell size
    AOA_threshold = 3*PDWParams.AOA.resolution;  % Threshold in degrees for AOA similarity
    Freq_threshold = 3*PDWParams.Freq.resolution;  % Threshold in Hz for frequency similarity

    for i = 1:length(PDWStruct)
        pdw = PDWStruct(i);
        added_to_cluster = false;
        
        % Check if this PDW fits into an existing cluster
        for c = 1:length(clusters)
            % Get the last PDW in the current cluster
            last_pdw = clusters{c}(end);
            
            % Check if the current PDW is similar enough to be added to this cluster
            if abs(pdw.AOA - last_pdw.AOA) <= AOA_threshold && ...
               abs(pdw.Freq - last_pdw.Freq) <= Freq_threshold
               
                % Add the PDW to this cluster
                clusters{c} = [clusters{c}, pdw];
                cell_size(c) = cell_size(c) + 1;  % Update the cluster size
                added_to_cluster = true;
                break;
            end
        end
        
        % If the PDW doesn't fit into any existing cluster, create a new one
        if ~added_to_cluster
            clusters{end+1} = pdw;
            cell_size(end+1) = 1;  % Initialize cell size for the new cluster
        end
    end
    
    if scatterPlot == 1
        figure;
        hold on;
        colororder("gem");
        for i = 1:length(clusters)
            frequency = [cell2mat(clusters(i)).Freq];
            pw = [cell2mat(clusters(i)).PW];
            aoa = [cell2mat(clusters(i)).AOA];
            scatter(frequency*1e-9,aoa);
            % Use PW plot to show why freq was used (less ambiguity)
        end
        hold off;
        xlabel('Frequency [GHz]');
        ylabel('AoA [degrees]');
        grid on;
        title('Scatter plot of clustered PDWs');
        % figureName = sprintf('/home/ryan/Documents/EEE5120Z/Project/Figures/Deint/clustering_scatter_plot.tex');
        % matlab2tikz('/home/ryan/Documents/EEE5120Z/Project/Figures/Deint/clustering_scatter_plot.tex');
        fprintf('Figure saved in Figures/ folder\n');
    end
    
    % Output the clusters to independent text files
    for c = 1:length(clusters)
        output_filename = sprintf('cluster_%d_PDW_output.txt', c);
        file_id = fopen(output_filename, 'w');
        
        % Write file header
        fprintf(file_id,'TOA,AOA,Amp,Freq,PW\n');

        for p = 1:length(clusters{c})
            pdw = clusters{c}(p);
            fprintf(file_id, '%.9f,%.5f,%d,%.5f,%.9f\n', ...
                    pdw.TOA, pdw.AOA, pdw.Amp, pdw.Freq, pdw.PW);
        end
        
        fclose(file_id);
        disp(['Cluster ', num2str(c), ' written to ', output_filename]);
    end
end