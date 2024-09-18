function plotAmplitudeTOA(filename,figureName,save)
    % Read the data from the CSV file
    data = readtable(filename);

    % Extract TOA and Amplitude columns
    TOA = data{:, 'TOA'};
    Amp_dBm = data{:, 'Amp'};

    % Create the plot
    figure;
    scatter(TOA, Amp_dBm);
    
    % Add labels and title
    xlabel('ToA [s]');
    ylabel('Amplitude [dBm]');
    title('Amplitude vs ToA');
    
    % Display grid for better visualization
    grid on;

    % If save is 1, save the figure in 'Figures/' directory
    if save == 1
        % Create the folder 'Figures/' if it doesn't exist
        if ~exist('Figures', 'dir')
            mkdir('Figures');
        end
        % cd('Figures/');
        
        % Save the figure as a high-quality PNG
        figureName = sprintf('/home/ryan/Documents/EEE5120Z/Project/Figures/PDW/%s.tex',figureName);
        matlab2tikz(figureName);
        fprintf('Figure saved as %s in Figures/ folder\n', figureName);
    end
end
