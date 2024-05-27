function [mdt, fr, tm] = trFNC_cluster_analysis(cluster_idx, cluster_num)
% TRFNC_ANALYSIS - Perform transition FNC (trFNC) analysis.
% Computes various metrics such as mean dwell time (mdt), fraction of
% recurrence (fr), and transition matrix (tm) for each subject. Additionally, it
% calculates group-level transition matrix (group_tm).
% Inputs:
%   idx: Matrix of cluster indices for each time point and subject.
%   cluster_num: Number of clusters.
% Outputs:
%   mdt: Mean dwell time for each subject and cluster.
%   fr: Fraction of recurrence for each subject and cluster.
%   tm: Transition matrix for each subject.

% Check input arguments
if nargin ~= 2
    error('Incorrect number of input arguments. Expected 2 inputs.');
end

if ~isnumeric(cluster_idx) || ~isnumeric(cluster_num)
    error('Inputs must be numeric.');
end

if size(cluster_idx, 2) < 2 || ~ismatrix(cluster_idx)
    error('idx must be a matrix with at least 2 columns.');
end

if ~isscalar(cluster_num) || cluster_num < 1 || mod(cluster_num, 1) ~= 0
    error('Invalid cluster_num. Must be a positive integer scalar.');
end

% Extract dimensions
subs = size(cluster_idx, 1);
tp = size(cluster_idx, 2);

% Initialize output variables
mdt = zeros(subs, cluster_num);
fr = zeros(subs, cluster_num);
tm = zeros(subs, cluster_num, cluster_num);

% Loop through subjects
for sub = 1:subs
    sub_idx = cluster_idx(sub, :);

    % Compute transition matrix
    temp_tm = compute_sub_tm(sub_idx, cluster_num);
    temp_tm(isnan(temp_tm)) = 0;
    tm(sub, :, :) = temp_tm;

    % Compute mean dwell time and fraction of recurrence
    for kn = 1:cluster_num
        mdt(sub, kn) = sum(sub_idx == kn) / compute_dwells(sub_idx, kn);
        fr(sub, kn) = sum(sub_idx == kn) / tp;
    end
end

end
