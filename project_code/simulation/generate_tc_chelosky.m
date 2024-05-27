function [Tc_tr, Tc_filt, Tc_uncorr, corr_tc, var_vec] = generate_tc_chelosky(options)
    % generate_tc_chelosky generates time series with specified correlation and variance characteristics
    %
    % Inputs:
    %   options - structure with the following fields:
    %       .band     - activity LPF passband [low high]
    %       .L        - time series length (number of points)
    %       .Tr       - time series sampling interval (in sec)
    %       .corr_tc  - vector of correlations at each time point
    %       .vartype  - 'equal' or 'random', specifies variance type
    %       .rngseed  - seed used to generate random values (0 means random)
    %
    % Outputs:
    %   Tc_tr    - time series with target correlation
    %   Tc_filt  - filtered time series without target correlation
    %   Tc_uncorr- unfiltered time series without target correlation
    %   corr_tc  - input correlation vector (returned for consistency)
    %   var_vec  - vector of variances at each time point
    %
    % Example usage:
    %   options.band = [0.01 0.15];
    %   options.L = 100;
    %   options.Tr = 1;
    %   options.corr_tc = ones(100, 1) * 0.5; % constant correlation of 0.5
    %   options.vartype = 'equal';
    %   options.rngseed = 0;
    %   [Tc_tr, Tc_filt, Tc_uncorr, corr_tc, var_vec] = generate_tc_chelosky(options);

    % Read the parameters
    band = options.band;
    L = options.L;
    Tr = options.Tr;
    corr_tc = options.corr_tc;
    vartype = options.vartype;
    rngseed = options.rngseed;
    
    % Derived parameters
    Fs = 1 / Tr;           % Sampling frequency
    nyquist = Fs / 2;      % Nyquist frequency
    Wp = band / nyquist;   % Normalized passband
    Ws = [band(1) * 0.5, band(2) * 2] / nyquist; % Stopband
    Rp = 3;                % Passband ripple in dB
    Rs = 30;               % Stopband attenuation in dB
    
    % Design Butterworth filter
    if all(Wp < 1)
        [n, Wn] = buttord(Wp, Ws, Rp, Rs);
        [kern_b, kern_a] = butter(n, Wn, 'bandpass');
    end

    % Time vector
    t = 0:Tr:(L - Tr);

    % Generate random variance if needed
    rand_var = rand(length(t), 1) / 2 + 1;

    % Set random seed if specified
    if rngseed ~= 0
        rng(rngseed);
    else
        rng('shuffle');
    end

    % Generate random time series
    pad_size = 40;
    tc1 = randn(length(t) + pad_size, 1); tc1 = tc1 - mean(tc1);
    tc2 = randn(length(t) + pad_size, 1); tc2 = tc2 - mean(tc2);

    % Mirror padding to avoid edge effects
    mirror_n = round(length(tc1) * 0.1);
    if all(Wp < 1)
        % Filter with mirroring to avoid edge effects
        tc1_m = [flip(tc1(1:mirror_n)); tc1; flip(tc1(end-mirror_n+1:end))];
        tc2_m = [flip(tc2(1:mirror_n)); tc2; flip(tc2(end-mirror_n+1:end))];
        
        tc1_filt = filtfilt(kern_b, kern_a, tc1_m);
        tc2_filt = filtfilt(kern_b, kern_a, tc2_m);
        
        tc1_filt = tc1_filt(mirror_n+1:end-mirror_n);
        tc2_filt = tc2_filt(mirror_n+1:end-mirror_n);
        
        tc1_filt = tc1_filt(pad_size/2+1:end-pad_size/2);
        tc2_filt = tc2_filt(pad_size/2+1:end-pad_size/2);
    else
        tc1_filt = tc1;
        tc2_filt = tc2;
    end
    
    % Unfiltered time series without target correlation
    Tc_uncorr = [tc1 tc2];
    Tc_filt = [tc1_filt tc2_filt];
    var_vec = zeros(L, 1);
    Tc_tr = zeros(L, 2);
    
    % Apply target correlation and variance
    for tn = 1:length(t)
        if strcmp(vartype, 'equal')
            var_vec(tn) = 1;
        elseif strcmp(vartype, 'random')
            var_vec(tn) = rand_var(tn);
        else
            error('Unknown vartype. Use "equal" or "random".');
        end
        cov_val = corr_tc(tn) * var_vec(tn);
        sigma = [var_vec(tn) cov_val; cov_val var_vec(tn)];
        sigma_chol = chol(sigma, 'lower');
        Tc_tr(tn, :) = [tc1_filt(tn) tc2_filt(tn)] * sigma_chol';
    end
end

