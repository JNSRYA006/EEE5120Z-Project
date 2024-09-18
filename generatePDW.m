function PDW = generatePDW(PDWParams,TOA,AOA,Amp, Freq, PW)

    % Hardcode whether or not to add uncertainty
    add_uncertainty = 1;

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
    TOA_discrete = discretiseVal(TOA,PDWParams.TOA.min,PDWParams.TOA.resolution,TOA_range);
    AOA_discrete = discretiseVal(AOA,PDWParams.AOA.min,PDWParams.AOA.resolution,AOA_range);
    Amp_discrete = discretiseVal(Amp,PDWParams.Amp.min,PDWParams.Amp.resolution,Amp_range);
    Freq_discrete = discretiseVal(Freq,PDWParams.Freq.min,PDWParams.Freq.resolution,Freq_range);
    PW_discrete = discretiseVal(PW,PDWParams.PW.min,PDWParams.PW.resolution,PW_range);
    % 2.5 - Add uncertainityto PDW values
    if add_uncertainty == 1
        PDW.TOA = addUncertainty(TOA_discrete,PDWParams.TOA.resolution);
        PDW.AOA = addUncertainty(AOA_discrete,PDWParams.AOA.resolution);
        PDW.Amp = addUncertainty(Amp_discrete,PDWParams.Amp.resolution);
        PDW.Freq = addUncertainty(Freq_discrete,PDWParams.Freq.resolution);
        PDW.PW = addUncertainty(PW_discrete,PDWParams.PW.resolution);
    else
        PDW.TOA = TOA_discrete;
        PDW.AOA = AOA_discrete;
        PDW.Amp = Amp_discrete;
        PDW.Freq = Freq_discrete;
        PDW.PW = PW_discrete;
    end
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

function uncertainValue = addUncertainty(value,resolution)
    chance = rand();
    if chance <= 0.1
        uncertainValue = value + resolution;
        % fprintf('Added uncertainity of %.4f to %.4f. New value = %.4f.\n',resolution,value,uncertainValue);
    elseif chance >= 0.9
        uncertainValue = value - resolution;
        % fprintf('Removed uncertainity of %.4f to %.4f. New value = %.4f.\n',resolution,value,uncertainValue);
    else
        uncertainValue = value;
    end
end