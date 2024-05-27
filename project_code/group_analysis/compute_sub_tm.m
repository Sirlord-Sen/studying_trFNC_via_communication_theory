function transition_prob = compute_sub_tm(sub_idx, kn_num)
% COMPUTE_SUB_FRACTION - Calculate transition probabilities for a given subject.
% Computes transition probabilities between clusters for a given subject,
% based on the cluster indices (sub_idx).
% Inputs:
%   sub_idx: Index vector representing cluster assignments over time for a subject.
%   kn_num: Number of clusters.
% Outputs:
%   transition_prob: Transition probabilities between clusters.

% Check input arguments
if nargin ~= 2
    error('Incorrect number of input arguments. Expected 2 inputs.');
end

if ~isnumeric(sub_idx) || ~isnumeric(kn_num)
    error('Inputs must be numeric.');
end

if ~isvector(sub_idx)
    error('sub_idx must be a vector.');
end

if ~isscalar(kn_num) || kn_num < 1 || mod(kn_num, 1) ~= 0
    error('Invalid kn_num. Must be a positive integer scalar.');
end

% Initialize transition counts and probabilities matrices
transition_cnts = zeros(kn_num, kn_num);
transition_prob = zeros(kn_num, kn_num);

% Count transitions
for i = 1:(length(sub_idx) - 1)
    from_state = sub_idx(i);
    to_state = sub_idx(i + 1);
    transition_cnts(from_state, to_state) = transition_cnts(from_state, to_state) + 1;
end

% Compute transition probabilities
for from_state = 1:kn_num
    ST = sum(sub_idx == from_state);  % Total occurrences of from_state
    for to_state = 1:kn_num
        transition_prob(from_state, to_state) = transition_cnts(from_state, to_state) / ST;
    end
end

end
