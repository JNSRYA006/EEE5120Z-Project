function PDWToCSV(PDW,filename)
    % Open file to write to
    file_id = fopen(filename,'w');

    % Check opening success
    if file_id == -1
        error('Could not open %s for writing.',filename);
    end
    
    % Write file header
    fprintf(file_id,'TOA,AOA,Amp,Freq,PW\n');

    % Write PDW to CSV file
    for i = 1:length(PDW)
        fprintf(file_id,'%.7f,%.5f,%d,%d,%.7f\n',PDW(i).TOA,PDW(i).AOA,PDW(i).Amp,PDW(i).Freq,PDW(i).PW);
    end

    fclose(file_id);
    fprintf('PDW written successfully to %s.\n',filename);
end