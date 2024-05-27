function dwell_count = compute_dwells(cluster_idx, kn_num)
% COMPUTE_DWELLS - Calculate the number of dwells of a specified cluster in a time series.
% Computes the number of dwells (consecutive occurrences) of a specified
% cluster (kn_num) in the given index vector (idx).
% Inputs:
%   idx: Index vector representing cluster assignments over time.
%   kn_num: Cluster number for which dwells are to be calculated.
% Outputs:
%   dwell_count: Number of dwells of the specified cluster.

% Check input arguments
if nargin ~= 2
    error('Incorrect number of input arguments. Expected 2 inputs.');
end

if ~isnumeric(cluster_idx) || ~isnumeric(kn_num)
    error('Inputs must be numeric.');
end

if ~isvector(cluster_idx)
    error('idx must be a vector.');
end

if ~isscalar(kn_num) || kn_num < 1 || mod(kn_num, 1) ~= 0
    error('Invalid kn_num. Must be a positive integer scalar.');
end

% Initialize variables
count = 0;      % Counter for consecutive occurrences
dwell_count = 0;    % Counter for the number of dwells

% Iterate through the vector
for i = 1:length(cluster_idx)
    if cluster_idx(i) == kn_num
        count = count + 1;
    else
        % If the value is not kn_num, reset the count
        count = 0;
    end
    
    % Check if we have found a complete dwell of kn_num
    if count == 1
        dwell_count = dwell_count + 1;
    end
end

end
