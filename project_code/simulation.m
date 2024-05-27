clc; clear; close all;

% Parameters
Tr = 1; % Sampling Time
Fs = 1/Tr; % Sampling frequency
nyquist = Fs/2; % Nyquist frequency
band = [0.01 0.15]; % Filtering frequency band
ph_band = [0.01 0.15]; % Filtering frequency band of phase
amp_band = [0.01 0.03]; % Filtering frequency band of amplitude
L = 1200; % Time length
t = (0:1/Fs:(Tr*(L)-1/Fs))'; % Time vector

% Scenarios
cases = ["crp_wins", "swpc_wins", "both_win"];

% Non-linear Effects
frequency = 0.1; % Frequency of the signals
pad_size = 100;

%Generate time-resolved connectivity for all scenarios
max_amp = [0.3, 0.9, 0.6];
all_connectivity_bands = [[0.02 0.03]; [0.002 0.004]; [0.005 0.008]];
all_connectivity_tcs = generate_connectivity_tcs(Tr, L, pad_size, max_amp, all_connectivity_bands);

% Options for Cholesky decomposition
options.L = L;                % time series length (in sec)
options.Tr = Tr;              % time series sampling freq
options.rngseed = 0;          % seed used to generate conn (0 means random)
options.vartype = 'equal';    % variance value at each time point.

% Simulation parameters
samps = 100; % Number of random samples
all_window_sizes = [15, 45, 75, 105, 135]';
all_bandwidths = [[0.01 0.04]; [0.03 0.07]; [0.05 0.09]; [0.07 0.11]; [0.09 0.13]];

% Initialization of metrics
swpc_corr = zeros(samps, length(all_window_sizes));
swpc_rmse = zeros(samps, length(all_window_sizes));
swpc_all = zeros(samps, length(all_window_sizes), L);
crp_corr = zeros(samps, length(all_bandwidths));
crp_rmse = zeros(samps, length(all_bandwidths));
crp_all = zeros(samps, length(all_bandwidths), L);

% Initialization of metrics for different cases
swpc_corr_all = zeros(samps, length(all_window_sizes), length(cases));
crp_corr_all = zeros(samps, length(all_bandwidths), length(cases));
swpc_rmse_all = zeros(samps, length(all_window_sizes), length(cases));
crp_rmse_all = zeros(samps, length(all_bandwidths), length(cases));

% Noise level
dbSNR = 50; % Signal to noise ratio in decibels

for samp = 1 : samps
    fprintf('Sample no: %d\n', samp)
    
    for caseIdx = 1:length(cases)
        scenario = cases(caseIdx);
        switch scenario
            case 'crp_wins'
                % Case 1: Ground truth in phase:Generate random amplitudes.
                connectivity_tc = all_connectivity_tcs(:, 1);
                
                tc1_amplitude = generate_random_ts([0.01 0.03], Tr, L);
                tc2_amplitude = generate_random_ts([0.01 0.03], Tr, L);
                tc1_amplitude = zscore(tc1_amplitude);
                tc2_amplitude = zscore(tc2_amplitude);
    
                % Generate reference phase
                ref_phase = generate_random_phase_ts([0.005 0.008], Tr, L, 1*pi/1);
                
                % Generate signals with phase and amplitude modulation
                tc1 = tc1_amplitude + 1*(sin(2 * pi * frequency * t + ref_phase + connectivity_tc));
                tc2 = tc2_amplitude + 1*(cos(2 * pi * frequency * t + ref_phase));

                % Set title
                title_text = 'CRP > SWPC';
                dbSNR = 50; %Signal to noise ratio in decibels
    
            case 'swpc_wins' 
                %Case 2: Generate full time series using Cholesky
                connectivity_tc = all_connectivity_tcs(:, 2);                
                cholesky_band = [0.01 0.15];

                options.corr_tc = connectivity_tc;    % Connectivity freq
                options.band = cholesky_band; % frequency band of signals in hertz

                [Tc_tr] = generate_tc_chelosky(options);
                tc1 = Tc_tr(:, 1);
                tc2 = Tc_tr(:, 2);

                % Set title
                title_text = 'SWPC > CRP';
                dbSNR = 10; %Signal to noise ratio in decibels

            case 'both_win'
                % Case 3: Use generate_tc_chelosky to generate amplitudes
                connectivity_tc = all_connectivity_tcs(:, 3);
                cholesky_band = [0.01 0.03];

                options.corr_tc = connectivity_tc;    % Connectivity freq
                options.band = cholesky_band; % frequency band of signals in hertz
    
                [Tc_tr] = generate_tc_chelosky(options);
                tc1_amplitude = Tc_tr(:, 1);
                tc2_amplitude = Tc_tr(:, 2);
    
                % Generate reference phase
                ref_phase = generate_random_phase_ts([0.005 0.008], Tr, L, 1*pi/1);
                
                % Generate signals with phase and amplitude modulation
                tc1 = tc1_amplitude + sin(2 * pi * frequency * t + ref_phase + connectivity_tc);
                tc2 = tc2_amplitude + cos(2 * pi * frequency * t + ref_phase);

                dbSNR = 10; %Signal to noise ratio in decibels

                % Set title
                title_text = 'SWPC ~= CRP';
        end

        %Add noise to signals
        tc1 = add_noise_to_signal(tc1, dbSNR);
        tc2 = add_noise_to_signal(tc2, dbSNR);

        for ws = 1:length(all_window_sizes)
            win_size = all_window_sizes(ws);
            
            % SWPC cropping snippet
            swpc_start = median(1:win_size);
            swpc_crop = (swpc_start:L-swpc_start-1)';
            
            % Compute sliding window Pearson correlation analysis
            swpc_tmp = sliding_window_pearson_correlation(L, win_size, tc1, tc2);
            swpc_connectivity_tc = connectivity_tc;
            swpc_all(samp, ws, :) = swpc_tmp;
            swpc_corr(samp, ws) = corr2(swpc_tmp(swpc_crop), swpc_connectivity_tc(swpc_crop));
            swpc_corr_all(samp, ws, caseIdx) = swpc_corr(samp, ws);
            swpc_rmse(samp, ws) = sqrt(mse(swpc_tmp(swpc_crop), swpc_connectivity_tc(swpc_crop)));
            swpc_rmse_all(samp, ws, caseIdx) = swpc_rmse(samp, ws);
        end
        
        for bw = 1:length(all_bandwidths)
            % Filtering
            tc1_h = bandpass_filtering(Tr, all_bandwidths(bw, :), tc1, pad_size);
            tc2_h = bandpass_filtering(Tr, all_bandwidths(bw, :), tc2, pad_size);
            
            % Compute cosine of relative phase
            crp_tmp = phase_synchrony(tc1_h, tc2_h);
            crp_all(samp, bw, :) = crp_tmp;
            crp_corr(samp, bw) = corr2(crp_tmp, connectivity_tc);
            crp_corr_all(samp, bw, caseIdx) = crp_corr(samp, bw);
            crp_rmse(samp, bw) = sqrt(mse(crp_tmp, connectivity_tc));
            crp_rmse_all(samp, bw, caseIdx) = crp_rmse(samp, bw);
        end
        
        % Plot connectivity_tc
        subplot(length(cases), 3, (caseIdx - 1) * 3 + 1);
        plot(t, connectivity_tc, 'k', LineWidth=2);
        title('Ground truth connectivity');
        ylim([-1 1]);
        
        % Plot SWPC and CRP metrics
        subplot(length(cases), 3, (caseIdx - 1) * 3 + 2);
        colors = ['b', 'b', 'b', 'b', 'b', 'r', 'r', 'r', 'r', 'r'];
        labels = {'WS (15s)','WS (45s)','WS (75s)','WS (105s)','WS (135s)'...
            ,'BW (0.01-0.04hz)', 'BW (0.03-0.07hz)', 'BW (0.05-0.9hz)', ...
            'BW (0.07-0.11hz)', 'BW (0.09-0.13hz)'};
        boxplot([swpc_corr_all(:, :, caseIdx), crp_corr_all(:, :, caseIdx)], 'Symbol', '+k' ,'Colors', colors ,'Labels', labels);
        title(['Correlation ' title_text]);
        ylim([-0.2 1]);

        % Plot SWPC and CRP metrics
        subplot(length(cases), 3, (caseIdx - 1) * 3 + 3);
        colors = ['b', 'b', 'b', 'b', 'b', 'r', 'r', 'r', 'r', 'r'];
        labels = {'WS (15s)','WS (45s)','WS (75s)','WS (105s)','WS (135s)'...
            ,'BW (0.01-0.04hz)', 'BW (0.03-0.07hz)', 'BW (0.05-0.9hz)', ...
            'BW (0.07-0.11hz)', 'BW (0.09-0.13hz)'};
        boxplot([swpc_rmse_all(:, :, caseIdx), crp_rmse_all(:, :, caseIdx)], 'Symbol', '+k' ,'Colors', colors ,'Labels', labels);
        title(['RMSE ' title_text]);
        ylim([0 1]);
    end
    
end

