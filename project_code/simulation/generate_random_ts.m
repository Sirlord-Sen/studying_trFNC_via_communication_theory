function ts_random = generate_random_ts(band, Tr, L)
    % generate_random_ts generates a bandpass-filtered random time series
    %
    % Inputs:
    %   band - bandpass filter frequencies [low high]
    %   Tr   - time series sampling interval (in sec)
    %   L    - length of the time series
    %
    % Output:
    %   ts_random - bandpass-filtered random time series with all positive values
    
    % Generate a random time series of length L
    ts = rand(L, 1);
    
    % Apply bandpass filtering to the time series
    ts_filtered = bandpass_filtering(Tr, band, ts);
    
    % Shift the filtered time series to ensure all values are positive
    ts_random = abs(min(ts_filtered)) + ts_filtered;
end