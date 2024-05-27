function [correlation_bins, subTcs_bin_idx] = swpc_ps_correlation_bins(swpc_ps_rho, bin_num)
    % swpc_ps_correlation_bins divides SWPC-PS correlation values into specified number of bins.
    %
    % Inputs:
    %   swpc_ps_rho:    2D array containing SWPC-PS correlation values (Subjects x Timepoints)
    %   bin_num:        Number of bins to divide the correlation values
    %
    % Outputs:
    %   correlation_bins:   2D array containing the bin ranges for correlation values
    %   subTcs_bin_idx:     3D array indicating the bin index for each subject's SWPC-PS correlation values

    % Get dimensions
    [subjects, timepoints] = size(swpc_ps_rho);

    % Pre-allocate arrays
    subTcs_bin_idx = zeros(bin_num, subjects, timepoints); % SWPC bin indices for each subject
    correlation_bins = zeros(bin_num, 2); % Bins for correlation values

    % Loop through different percentile ranges
    for p = 1:bin_num
        fprintf(['bin: ' num2str(p) '/' num2str(bin_num) '\n']);
        
        % Calculate bin range
        lower_bin = (p - 1) / bin_num;
        upper_bin = p / bin_num;
        corr_bin = [lower_bin upper_bin];
        correlation_bins(p, :) = corr_bin;

        % Assign SWPC-PS correlation values to bins for each subject
        for sub = 1:subjects
            swpc_ps_sub = swpc_ps_rho(sub, :);
            sub_indices = (swpc_ps_sub >= corr_bin(1) & swpc_ps_sub <= corr_bin(2));
            subTcs_bin_idx(p, sub, :) = sub_indices;
        end
    end
end
