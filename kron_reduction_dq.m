function Y_reduced = kron_reduction_dq(Y, eliminate_buses)
%KRON_REDUCTION_DQ Performs Kron reduction for dq-frame admittance matrices
%   Y_reduced = kron_reduction_dq(Y, eliminate_buses)
%   
%   Inputs:
%   Y - N x N admittance matrix in dq coordinates (N must be even)
%   eliminate_buses - List of bus numbers (1-based) to eliminate
%
%   Output:
%   Y_reduced - Reduced admittance matrix keeping dq components of remaining buses

% Validate matrix dimensions
[N, M] = size(Y);
if N ~= M || mod(N, 2) ~= 0
    error('Y must be a square matrix with even dimensions (2n x 2n)');
end
n_buses = N/2;

% Validate bus numbers
eliminate_buses = unique(eliminate_buses);
if any(eliminate_buses < 1) || any(eliminate_buses > n_buses)
    error('Invalid bus numbers specified');
end

% Convert bus numbers to node indices (2 nodes per bus: d and q)
eliminate_nodes = [];
for bus = eliminate_buses
    base_idx = (bus-1)*2;
    eliminate_nodes = [eliminate_nodes, base_idx+1, base_idx+2];
end

% Get retained nodes (automatically keeps dq pairs together)
all_nodes = 1:N;
retained_nodes = setdiff(all_nodes, eliminate_nodes);

% Reorder matrix: retained nodes first, eliminated last
new_order = [retained_nodes, eliminate_nodes];
Y_ordered = Y(new_order, new_order);

% Partition the matrix
n_retained = length(retained_nodes);
Y_rr = Y_ordered(1:n_retained, 1:n_retained);%  Y_rr = minreal(Y_rr);
Y_re = Y_ordered(1:n_retained, n_retained+1:end);% Y_re = minreal(Y_re);
Y_er = Y_ordered(n_retained+1:end, 1:n_retained);% Y_er = minreal(Y_er);
Y_ee = Y_ordered(n_retained+1:end, n_retained+1:end);% Y_ee = minreal(Y_ee);
 
% Check invertibility and compute Schur complement
try
    Y_ee_inv = inv(Y_ee);
catch
    error('Y_ee matrix is singular - Kron reduction not possible');
end

% Calculate reduced admittance matrix
Y_reduced = Y_rr - Y_re * Y_ee_inv * Y_er;

% Maintain original dq-pair ordering in output
[~, sort_idx] = sort(retained_nodes);
Y_reduced = Y_reduced(sort_idx, sort_idx);

end