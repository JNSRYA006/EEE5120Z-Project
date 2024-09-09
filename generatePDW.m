function PDW = generatePDW(PDWParams,TOA,AOA,Amp, Freq, PW)

    % Check if input values are in the allowable range
    try
        % TOA
        if TOA < PDWParams.TOA.min || TOA > PDWParams.TOA.max
            error('TOA must be within the range %.0f to %.2f seconds',PDWParams.TOA.min,PDWParams.TOA.max);
        end
    catch ME
        error('TOA Error: %s\n', ME.message);
    end
    
    try
        % AOA
        if AOA < PDWParams.AOA.min || AOA > PDWParams.AOA.max
            error('AOA must be within the range %.0f to %.2f degrees',PDWParams.AOA.min,PDWParams.AOA.max);
        end
    catch ME
        error('AOA Error: %s\n', ME.message);
    end
    
    try
        % Amplitude
        % disp(Amp);
        if Amp < PDWParams.Amp.min || Amp > PDWParams.Amp.max
            error('Amplitude must be within the range %.2f to %.2f dB',PDWParams.Amp.min,PDWParams.Amp.max);
        end
    catch ME
        error('Amplitude Error: %s\n', ME.message);
    end
    
    try
        % Frequency
        if Freq < PDWParams.Freq.min || Freq > PDWParams.Freq.max
            error('Frequency must be within the range %.2f to %.2f GHz',PDWParams.Freq.min*1e-9,PDWParams.Freq.max*1e-9);
        end
    catch ME
        error('Frequency Error: %s\n', ME.message);
    end
    
    try
        % Pulse Width
        if PW < PDWParams.PW.min || PW > PDWParams.PW.max
            error('Pulse Width must be within the range %.2f to %.2f µs',PDWParams.PW.min*1e6, PDWParams.PW.max*1e6);
        end
    catch ME
        error('Pulse Width Error: %s\n', ME.message);
    end

    % Create the PDW structure by discretising the input parameters
    % Will be the closet value
    % Calculate range of values
    TOA_range = calcRange(PDWParams.TOA.bits,PDWParams.TOA.resolution);
    AOA_range = calcRange(PDWParams.AOA.bits,PDWParams.AOA.resolution);
    Amp_range = calcRange(PDWParams.Amp.bits,PDWParams.Amp.resolution);
    Freq_range = calcRange(PDWParams.Freq.bits,PDWParams.Freq.resolution);
    PW_range = calcRange(PDWParams.PW.bits,PDWParams.PW.resolution);
    % Discretise
    PDW.TOA = discretiseVal(TOA,PDWParams.TOA.min,PDWParams.TOA.resolution,TOA_range);
    PDW.AOA = discretiseVal(AOA,PDWParams.AOA.min,PDWParams.AOA.resolution,AOA_range);
    PDW.Amp = discretiseVal(Amp,PDWParams.Amp.min,PDWParams.Amp.resolution,Amp_range);
    PDW.Freq = discretiseVal(Freq,PDWParams.Freq.min,PDWParams.Freq.resolution,Freq_range);
    PDW.PW = discretiseVal(PW,PDWParams.PW.min,PDWParams.PW.resolution,PW_range);
end

function discretised = discretiseVal(value,min_val,resolution,range)
    norm_val = (value-min_val)/range;
    levels = range/resolution;
    discrete_level = round(norm_val*levels);
    discretised = min_val + discrete_level*resolution;
    % range = (2^bits -1)*resolution;
end

function range = calcRange(bits,resolution)
    range = (2^bits -1)*resolution;
end