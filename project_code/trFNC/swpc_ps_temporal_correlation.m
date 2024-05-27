function [swpc_ps_rho] = swpc_ps_temporal_correlation(swpc_est, ps_est)
    % swpc_ps_temporal_correlation calculates temporal correlation between PS estimates and SWPC estimates.
    %
    % Inputs:
    %   swpc_est:   4D array containing SWPC estimates (Subjects x Time x Regions x Regions)
    %   ps_est:     4D array containing PS estimates (Subjects x Time x Regions x Frequencies)
    %
    % Output:
    %   swpc_ps_rho:    2D array containing temporal correlation between PS and SWPC estimates for each subject

    % Pre-allocate arrays
    swpc_ps_rho = zeros(size(ps_est, 1), size(ps_est, 2)); % Concatenated time series per subject
    rho_no_padding = zeros(size(ps_est, 1), size(ps_est, 2), size(ps_est, 2)); % 3D array without padded zeros

    % Loop over subjects
    for sub = 1:size(ps_est, 1)
        fprintf('sub: %d\n', sub);
        sub_iPS_est = squeeze(ps_est(sub, :, :, :)); % PS estimate converted to 3D
        sub_SWPC_est = squeeze(swpc_est(sub, :, :, :)); % SWPC estimate converted to 3D
        
        % Compute Spearman correlation without padded zeros
        rho_no_padding(sub, :, :) = corr(icatb_mat2vec(sub_iPS_est)', icatb_mat2vec(sub_SWPC_est)', 'Type', 'Spearman');

        % Extract diagonal values
        swpc_ps_rho(sub, :) = diag(squeeze(rho_no_padding(sub, :, :)));
    end
end
