function phase_tc = generate_random_phase_ts(band, Tr, L, max_amp)
    % generate_random_phase_ts generates a time series with random phase within a specified frequency band
    %
    % Inputs:
    %   band    - bandpass filter frequencies [low high]
    %   Tr      - time series sampling interval (in sec)
    %   L       - length of the time series
    %   max_amp - maximum amplitude for the time series
    %
    % Output:
    %   phase_tc - generated time series with random phase and specified maximum amplitude
    
    % Generate random initial time series of half the desired length
    ts = rand(L / 2 + 1, 1);
    
    % Despike the time series
    ts = icatb_despike_tc(ts, Tr);
    
    % Apply bandpass filtering to the time series
    ts_filtered = bandpass_filtering(Tr, band, ts);
    
    % Fix the amplitude
    amp = max(abs(ts_filtered));
    ts_filtered = ts_filtered * (max_amp / amp);
    
    % Make symmetric to ensure the signal is real
    phase_tc = [ts_filtered' (ts_filtered(L / 2:-1:2))']';
end