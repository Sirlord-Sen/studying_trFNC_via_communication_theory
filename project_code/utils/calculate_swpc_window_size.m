function win_size = calculate_swpc_window_size(Tr, band)
    % calculate_swpc_window_size computes the optimal window size for SWPC.
    %
    % Inputs:
    %   Tr:     Repetition Time (sampling interval)
    %   band:   Frequency band for filtering
    %
    % Outputs:
    %   win_size:   Optimal window size

    win_size = sqrt(((0.88 / band(1)) * (1 / Tr))^2 + 1);
    win_size = ceil(win_size); % Round number up
    % Ensure window size is odd
    if mod(win_size, 2) == 0
        win_size = win_size + 1;
    end
end
