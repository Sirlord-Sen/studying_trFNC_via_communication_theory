function [p_values, fdr_p_values, s_values, sign_values] = compute_cluster_stats(sub_centroids_group1, sub_centroids_group2, cluster_num)
    % compute_cluster_stats Computes statistical metrics for cluster data
    %
    % This function calculates p-values, FDR p-values, s-values, and sign values for given
    % group1 and group2 sub-centroid data over a specified number of clusters.
    %
    % INPUTS:
    %   sub_centroids_group1: A 3D numeric array of sub-centroids for group 1 (subjects x clusters x features).
    %   sub_centroids_group2: A 3D numeric array of sub-centroids for group 2 (subjects x clusters x features).
    %   cluster_num: An integer specifying the number of clusters.
    %
    % OUTPUTS:
    %   p_values: A matrix of p-values for each cluster and feature.
    %   fdr_p_values: A matrix of FDR q-values (corrected p-values) for each cluster and feature.
    %   s_values: A matrix of s-values (signed -log10 of FDR p-values) for each cluster and feature.
    %   sign_values: A matrix of sign values indicating the direction of the mean difference.
    %
    % Example:
    %   [p_values, fdr_q_values, s_values, sign_values] = compute_cluster_stats(sub_centroids_group1, sub_centroids_group2, 5);

    % Input validation
    if nargin ~= 3
        error('Incorrect number of inputs. Expected 3 inputs: sub_centroids_group1, sub_centroids_group2, and cluster_num.');
    end
    
    if ~isnumeric(sub_centroids_group1) || ~isnumeric(sub_centroids_group2) || ~isnumeric(cluster_num)
        error('All inputs must be numeric.');
    end
    
    if ndims(sub_centroids_group1) ~= 3 || ndims(sub_centroids_group2) ~= 3
        error('sub_centroids_group1 and sub_centroids_group2 must be 3-dimensional arrays.');
    end
    
    if size(sub_centroids_group1, 2) ~= cluster_num || size(sub_centroids_group2, 2) ~= cluster_num
        error('The second dimension of sub_centroids_group1 and sub_centroids_group2 must match cluster_num.');
    end
    
    % Initialize variables
    feat_num = size(sub_centroids_group1, 3);
    p_values = zeros(cluster_num, feat_num);
    fdr_p_values = zeros(cluster_num, feat_num);
    s_values = zeros(cluster_num, feat_num);
    sign_values = zeros(cluster_num, feat_num);
    
    % Computation loop
    for kn = 1:cluster_num
        for k = 1:feat_num
            group1_vec = sub_centroids_group1(:, kn, k);
            group2_vec = sub_centroids_group2(:, kn, k);
            
            group1_feat = group1_vec(~isnan(group1_vec));
            group2_feat = group2_vec(~isnan(group2_vec));
            
            sign_values(kn, k) = sign(mean(group1_feat) - mean(group2_feat));
            
            if isempty(group2_feat) || isempty(group1_feat)
                p_values(kn, k) = NaN;
            else
                [~, p_values(kn, k)] = ttest2(group2_feat, group1_feat);
            end
        end
        
        fdr_p_values(kn, :) = mafdr(p_values(kn, :), 'Bhfdr', true);
        s_values(kn, :) = -log10(fdr_p_values(kn, :)) .* sign_values(kn, :);
    end
    
    % Convert to matrix format
    p_values = icatb_vec2mat(p_values);
    fdr_p_values = icatb_vec2mat(fdr_p_values);
    s_values = icatb_vec2mat(s_values);
    sign_values = icatb_vec2mat(sign_values);
end
