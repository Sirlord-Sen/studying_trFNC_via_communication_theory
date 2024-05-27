%Add utils path
clc;clear;close all
addpath("./utils/");
addpath("./trFNC/");
%% Data
%Example of post-processed fMRI timecourses
subTcs = randn(5, 150, 10);
Tr = 2; %sampling time in seconds

%% SWPC postprocessing
swpc_band = [0.01 0.15]; % fMRI frequency band for SWPC: [low high]
cutoff_lim = [swpc_band(1)*0.7 swpc_band(2)*1.3];
display_flag = true;
subTcs_swpc = post_processing_subject_timecourses(subTcs, Tr, swpc_band, cutoff_lim, display_flag);

%% SWPC computation
%set window size to [] if you want window size to be calculated based on
%the -3dB point of the high pass filter. Otherwise, specify an odd window
%size
[swpc_zeros, swpc, win_center, win_size] = calculate_swpc(subTcs_swpc, Tr, swpc_band, [], "rectangular");

%% PS postprocessing
ps_band = [0.03 0.07]; % fMRI frequency band for PS: [low high]
cutoff_lim = [ps_band(1)*0.7 ps_band(2)*1.3];
display_flag = true;
subTcs_ps = post_processing_subject_timecourses(subTcs, Tr, ps_band, cutoff_lim, display_flag);

%% PS computation (No matching with SWPC)
[ps, ~] = calculate_ps(subTcs_ps);
%% K-means
cluster_num = 4;
kmeans_dist = 'cityblock'; 
rep = 20;
swpc_vectorized = icatb_mat2vec(reshape(swpc, size(swpc, 1) * size(swpc, 2), size(swpc, 3), size(swpc, 3)));
ps_vectorized = icatb_mat2vec(reshape(ps, size(ps, 1) * size(ps, 2), size(ps, 3), size(ps, 3)));

% SWPC clustering
[swpc_idx, swpc_c, swpc_sumd, swpc_D] = calculate_kmeans(swpc_vectorized, cluster_num, kmeans_dist, rep);

% PS clustering
[ps_idx, ps_c, ps_sumd, ps_D] = calculate_kmeans(ps_vectorized, cluster_num, kmeans_dist, rep);
