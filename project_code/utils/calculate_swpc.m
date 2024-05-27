function [swpc_est_zeros, swpc_est, win_center, win_size] = calculate_swpc(subTcs, Tr, band, win_size, window_type)
    % calculate_swpc computes Sliding Window Pearson Correlation (SWPC) matrices.
    %
    % Inputs:
    %   subTcs:         3D array containing time series data for all subjects (Subjects x Time x Regions)
    %   Tr:             Repetition Time (sampling interval)
    %   band:           Frequency band for filtering
    %   win_size:       Size of the sliding window (optional, if not provided, it will be calculated)
    %   window_type:    Type of windowing function ('rectangular', 'hamming', 'gaussian')
    %
    % Outputs:
    %   swpc_est_zeros: 4D array containing SWPC matrices with zero padding (Subjects x Time x Regions x Regions)
    %   swpc_est:       4D array containing SWPC matrices aligned to the center of the windows (Subjects x Time x Regions x Regions)
    %   win_center:     Vector containing the center index of each window
    %   win_size:       Size of the sliding window

    % If window size was not provided, compute it
    if nargin < 4 || isempty(win_size)
        win_size = calculate_swpc_window_size(Tr, band); % Get window size
    end

    % Ensure window size is odd
    if mod(win_size, 2) == 0
        win_size = win_size + 1;
        warning(['Window size has been set to: ', num2str(win_size)]);
    end

    % Adjust window size if it's less than 11
    if win_size < 11
        win_size = 11;
    end

    win_center = zeros(1, size(subTcs, 2) - win_size);
    swpc_est_zeros = zeros(size(subTcs, 1), size(subTcs, 2), size(subTcs, 3), size(subTcs, 3)); % Pre-allocate 4D SWPC matrix

    % Compute window based on the type provided
    switch window_type
        case 'rectangular'
            window = ones(win_size, 1);
        case 'hamming'
            window = hamming(win_size);
        case 'gaussian'
            window = gausswin(win_size);
        otherwise
            error('Invalid window type. Please choose from: rectangular, hamming, gaussian');
    end

    % Loop over subjects
    for sub = 1:size(subTcs, 1)
        fprintf("sub: %d\n", sub);
        widx = 1:win_size; % Define window 1-D matrix

        % Slide the window along the time dimension
        while widx(end) < size(subTcs, 2)
            % Apply window to subTcs data
            windowed_subTcs = bsxfun(@times, zscore(squeeze(subTcs(sub, widx, :))), window);
            % Compute Pearson correlation and store in SWPC_est_zeros
            swpc_est_zeros(sub, median(widx), :, :) = corr(windowed_subTcs);
            % Store center index of the window
            win_center(widx(1)) = median(widx);
            widx = widx + 1;
        end
    end

    % Align SWPC matrices to the center of the windows
    start = min(win_center);
    finish = max(win_center);
    swpc_est = swpc_est_zeros(:, start:finish, :, :);
end