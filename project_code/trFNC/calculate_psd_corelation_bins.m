function [subject_bin_psd, subject_bin_fft, subject_bin_subTc] = calculate_psd_corelation_bins(bin_num, subTcs_bin_idx, subTcs, Tr, win_size)
% CALCULATE_PRCTILE - Computes the power spectral density of the fMRI data corresponding to correlation bins.
% Inputs:
%   bin_num: Number of correlation bins.
%   subTcs_bin_idx: Indices of percentiles for each subject.
%   subTcs: fMRI time courses subject x time x nodes/components.
%   Tr: Sampling time.
%   win_size: Size of the window.
% Outputs:
%   subject_prctile_psd: Percentile Power Spectral Density data for each subject.
%   subject_prctile_fft: Percentile FFT data for each subject.
%   subject_prctile_subTc: Percentile subTc data for each subject.

% Check input arguments
if nargin ~= 5
    error('Incorrect number of input arguments. Expected 5 inputs.');
end

if ~isscalar(bin_num) || bin_num <= 0 || ~isnumeric(bin_num) || mod(bin_num, 1) ~= 0
    error('Invalid number of correlation bins. Must be a positive integer scalar.');
end

if ~isscalar(win_size) || win_size <= 0 || ~isnumeric(win_size) || mod(win_size, 1) ~= 0
    error('Invalid window size. Must be a positive integer scalar.');
end

if ~isscalar(Tr) || Tr <= 0 || ~isnumeric(Tr)
    error('Invalid sampling time. Must be a positive scalar.');
end

% Extract number of subjects and nodes/components
[subjects_cnt, ~, components_cnt] = size(subTcs);

%Sampling frequency
Fs = 1/Tr;

% Initialize output variables
subject_bin_psd = zeros(bin_num, subjects_cnt, ceil(win_size/2), components_cnt);
subject_bin_fft = zeros(bin_num, subjects_cnt, win_size, components_cnt);
subject_bin_subTc = zeros(bin_num, subjects_cnt, win_size, components_cnt);

% Loop through correlation bins
for p = 1:bin_num
    fprintf("bin: %d/%d \n", p, bin_num);
    % Loop through subjects
    for sub = 1:subjects_cnt
        % Get indices of non-zero values for current bin and subject
        ia = find(squeeze(subTcs_bin_idx(p, sub, :)));
        if ~isempty(ia)
            % Initialize temporary variables to store data for this subject
            % and bin
            tmp_fft = zeros(length(ia), win_size, components_cnt);
            tmp_psd = zeros(length(ia), ceil(win_size/2), components_cnt);
            tmp_subTc = zeros(length(ia), win_size, components_cnt);
            % Loop through indices and compute FFT and PSD
            for idx = 1:length(ia)
                % Get actual index for the window
                ia_actual = floor(win_size/2) + ia(idx);
                % Define the windowed indices
                windowed_idx = ia_actual - floor(win_size/2):ia_actual + floor(win_size/2);
                % Store the windowed subTc data
                tmp_subTc(idx, :, :) = squeeze(subTcs(sub, windowed_idx, :));
                % Compute FFT of the windowed subTc data
                tmp_fft(idx, :, :) = fft(squeeze(subTcs(sub, windowed_idx, :)));
                % Compute PSD of the windowed subTc data
                [tmp_psd(idx, :, :), ~] = periodogram(squeeze(subTcs(sub, windowed_idx, :)), [], win_size, Fs);
            end
            % Check if there is only one data point for this bin and subject
            if size(tmp_fft, 1) == 1
                % Store data directly
                subject_bin_subTc(p, sub, :, :) = squeeze(tmp_subTc);
                subject_bin_fft(p, sub, :, :) = squeeze(tmp_fft);
                subject_bin_psd(p, sub, :, :) = squeeze(tmp_psd);
            else
                % Compute mean and store data
                subject_bin_fft(p, sub, :, :) = squeeze(mean(tmp_fft));
                subject_bin_subTc(p, sub, :, :) = squeeze(mean(tmp_subTc));
                subject_bin_psd(p, sub, :, :) = squeeze(mean(tmp_psd));
            end
        else
            % If no data for this bin and subject, set PSD to NaN
            subject_bin_psd(p, sub, :, :) = nan(ceil(win_size/2), components_cnt);
        end
    end
end

end
