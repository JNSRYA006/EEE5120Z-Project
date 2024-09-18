function PDWStruct = importPDWLists(filename) 

    % Reads a CSV file containing PDW data and returns it as a struct array.
    %
    % Input:
    %   filename - The name or path of the CSV file to read.
    %
    % Output:
    %   pdw_struct - A struct array with fields TOA, AOA, Amp, Freq, and PW.
    
    % Read the CSV file into a table
    pdw_table = readtable(filename, 'Delimiter', ',', 'ReadVariableNames', true);
    
    % Convert the table to a struct array
    PDWStruct = table2struct(pdw_table, 'ToScalar', false)'; 
end