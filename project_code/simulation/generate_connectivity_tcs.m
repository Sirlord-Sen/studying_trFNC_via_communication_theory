function all_connectivity_tcs = generate_connectivity_tcs(Tr, L, pad_size, max_amp, all_connectivity_bands)
    % generate_connectivity_tcs generates ground truth connectivity time courses
    %
    % Inputs:
    %   Tr                    - time series sampling interval (in sec)
    %   L                     - length of the time series
    %   pad_size              - padding size for the filtering
    %   max_amp               - vector of maximum amplitudes for each connectivity band
    %   all_connectivity_bands - matrix of connectivity frequency bands (each row is a band)
    %
    % Output:
    %   all_connectivity_tcs - matrix of generated connectivity time courses
    
    % Initialize the output matrix
    num_cases = length(max_amp);
    all_connectivity_tcs = zeros(L, num_cases);
    
    % Generate ground truth connectivity time courses
    for cb = 1:num_cases
        temp_connectivity = bandpass_filtering(Tr, all_connectivity_bands(cb, :), rand(L, 1), pad_size);
        all_connectivity_tcs(:, cb) = temp_connectivity * (max_amp(cb) / max(abs(temp_connectivity)));
    end
end
