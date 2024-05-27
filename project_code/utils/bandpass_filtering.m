function tc = bandpass_filtering(Tr, band, tc, pad_size)
    % bandpass_filtering applies a Butterworth bandpass filter to a time series
    %
    % Inputs:
    %   Tr       - time series sampling interval (in sec)
    %   band     - bandpass filter frequencies [low high]
    %   tc       - input time series
    %   pad_size - number of points for padding to avoid edge effects
    %
    % Output:
    %   tc       - filtered time series

    % Set default pad_size if not provided
    if nargin < 4
        pad_size = 100;
    end
    
    % Sampling frequency
    Fs = 1 / Tr;
    nyquist = Fs / 2;
    
    % Normalized passband and stopband frequencies
    Wp = band / nyquist;
    Ws = [band(1) * 0.5, band(2) * 2] / nyquist;
    
    % Passband and stopband ripple/attenuation
    Rp = 3;  % Passband ripple in dB
    Rs = 30; % Stopband attenuation in dB
    
    % Design Butterworth bandpass filter
    if all(Wp < 1)
        [n, Wn] = buttord(Wp, Ws, Rp, Rs);
        [kern_b, kern_a] = butter(n, Wn, 'bandpass');
    end
    
    % Pad the time series with zeros
    tc = [zeros(pad_size / 2, 1); tc; zeros(pad_size / 2, 1)];

    % Mirror padding to avoid edge effects
    mirror_n = round(length(tc) * 0.1);
    if all(Wp < 1)
        % Create mirrored time series
        tc_m = [flip(tc(1:mirror_n)); tc; flip(tc(end-mirror_n+1:end))];
        
        % Apply the bandpass filter
        tc_filt = filtfilt(kern_b, kern_a, tc_m);
        
        % Remove the mirrored parts
        tc_filt = tc_filt(mirror_n + 1:end - mirror_n);
        
        % Remove the zero padding
        tc = tc_filt(pad_size / 2 + 1:end - pad_size / 2);
    end
end