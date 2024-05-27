function swpc = sliding_window_pearson_correlation(L, win_size, tc1, tc2)
    % sliding_window_pearson_correlation computes the sliding window Pearson correlation between two time series
    %
    % Inputs:
    %   L        - length of the output correlation vector
    %   win_size - size of the sliding window
    %   tc1      - first input time series
    %   tc2      - second input time series
    %
    % Output:
    %   swpc     - vector of sliding window Pearson correlations
    
    % Initialize the output vector
    swpc = zeros(L, 1);
    
    % Define the initial sliding window indices
    widx = 1:win_size;

    % Compute sliding window Pearson correlation
    while widx(end) <= length(tc1)
        % Calculate the correlation for the current window
        current_corr = corr(tc1(widx), tc2(widx));
        
        % Assign the correlation value to the middle of the current window
        swpc(median(widx)) = current_corr;
        
        % Shift the window by one time point
        widx = widx + 1;
    end
end