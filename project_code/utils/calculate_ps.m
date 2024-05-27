function [crp_est, srp_est] = calculate_ps(subTcs)
    % calculate_crp computes Cosine of Relative Phase (CRP) and Sine of Relative Phase (SRP) matrices.
    %
    % Inputs:
    %   subTcs:     3D array containing time series data for all subjects (Subjects x Time x Regions)
    %
    % Outputs:
    %   crp_est:    4D array containing CRP matrices (Subjects x Time x Regions x Regions)
    %   srp_est:    4D array containing SRP matrices (Subjects x Time x Regions x Regions)

    % Pre-allocate matrices
    crp_est = zeros(size(subTcs, 1), size(subTcs, 2), size(subTcs, 3), size(subTcs, 3));
    srp_est = zeros(size(subTcs, 1), size(subTcs, 2), size(subTcs, 3), size(subTcs, 3));

    % Loop over subjects
    for sub = 1:size(subTcs, 1)
        fprintf("sub: %d\n", sub);
        
        % Calculate instantaneous phase using Hilbert transform
        inst_phase = angle(hilbert(squeeze(subTcs(sub, :, :))));
        
        % Loop over time points
        for t = 1:size(crp_est, 2)
            % Calculate cosine and sine of relative phase matrices
            cos_ip = cos(bsxfun(@minus, inst_phase(t, :)', inst_phase(t, :)));
            srp_ip = sin(bsxfun(@minus, inst_phase(t, :)', inst_phase(t, :)));
            
            % Ensure there are no NaN or Inf values in matrices
            cos_ip(isnan(cos_ip)) = 0;
            cos_ip(isinf(cos_ip)) = 0;
            srp_ip(isnan(srp_ip)) = 0;
            srp_ip(isinf(srp_ip)) = 0;
            
            % Store CRP and SRP matrices
            crp_est(sub, t, :, :) = cos_ip;
            srp_est(sub, t, :, :) = srp_ip;
        end
    end
end
