%Add all paths
clc;clear;close all
addpath("./trFNC/");
addpath("./utils/");
%% Data
%Example of post-processed fMRI timecourses
subTcs = randn(5, 150, 10);
Tr = 2; %sampling time in seconds

%% Set case
case_num = 1;
if case_num == 1
    swpc_band = [0.01 0.15]; % fMRI frequency band for SWPC: [low high]
    ps_band = [0.03 0.07]; % fMRI frequency band for PS: [low high]
    win_size = 88/Tr;
elseif case_num == 2
    swpc_band = [0.03 0.07]; % fMRI frequency band for SWPC: [low high]
    ps_band = [0.03 0.07]; % fMRI frequency band for PS: [low high]
    win_size = 30/Tr;
elseif case_num == 3
    swpc_band = [0.03 0.07]; % fMRI frequency band for SWPC: [low high]
    ps_band = [0.03 0.07]; % fMRI frequency band for PS: [low high]
    win_size = 88/Tr;
end

%% SWPC postprocessing
cutoff_lim = [swpc_band(1)*0.7 swpc_band(2)*1.3];
display_flag = true;
subTcs_swpc = post_processing_subject_timecourses(subTcs, Tr, swpc_band, cutoff_lim, display_flag);

%% SWPC computation
%set window size to [] if you want window size to be calculated based on
%the -3dB point of the high pass filter. Otherwise, specify an odd window
%size
[swpc_zeros, swpc, win_center, win_size] = calculate_swpc(subTcs_swpc, Tr, swpc_band, win_size, "rectangular");

%% PS postprocessing
cutoff_lim = [ps_band(1)*0.7 ps_band(2)*1.3];
display_flag = true;
subTcs_ps = post_processing_subject_timecourses(subTcs, Tr, ps_band, cutoff_lim, display_flag);

%% PS computation
[ps, ~] = calculate_ps(subTcs_ps);
ps = ps(:, win_center, :, :); % Matching PS to SWPC
%% Spearman correlation between SWPC and PS
swpc_ps_rho = swpc_ps_temporal_correlation(swpc, ps);
histogram(swpc_ps_rho)
%% Splitting temporally matched indexes into correlation bins
bin_num = 5; % 5 correlation bins
[correlation_bins, subTcs_bin_idx] = swpc_ps_correlation_bins(swpc_ps_rho, bin_num);

%% Get subject PSD estimations per correlation bin
[subject_bin_psd, subject_bin_fft, subject_bin_subTc] = calculate_psd_corelation_bins(bin_num, subTcs_bin_idx, subTcs, Tr, win_size);
