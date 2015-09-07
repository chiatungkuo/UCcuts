function [V, P] = multisp(Gs, k, alpha)
%
% Implements the multiview spectral clustering proposed by Zhou et. al. in
% ICML 2007.
%
% Inputs
%   Gs: cell array containing graphs (affinity matrices) of the same sizes
%   k: number of clusters (i.e. number of eigenvectors to return)
%   alpha (optional): relative weights in combining graphs; vector of length 
%       the number of graphs
% Outputs
%   V: n by k matrix of columnwise stacked eigenvectors (hard clusters could be
%       determined upon user's choice.) 
%

nG = length(Gs);        % number of graphs
n = size(Gs{1}, 1);     % size of graphs
if nargin < 3, alpha = ones(1, nG)/nG; end

beta = zeros(n, nG);    
station = zeros(n, nG); % stationary distribution
Ps = cell(1, nG);       % transition probability matrices
ssum = zeros(n, 1);     % denominator in computing beta

for i = 1:nG
    A = Gs{i};
    d = sum(A, 2);      % degree vector
    Ps{i} = diag(1./d)*A;   
    vol = sum(d);       % volume of the graph
    station(:, i) = d/vol;      
    ssum = ssum + alpha(i)*station(:, i);
end

% compute beta and mixture probability transition
P = zeros(n, n);
for i = 1:nG
    beta(:, i) = alpha(i)*station(:, i)./ssum;
    P = P + diag(beta(:, i))*Ps{i};
end

% compute generalized eigenvectors, (real-valued) relaxed solution
Z = diag(ssum);
L = Z - (Z*P + P'*Z)/2;
[V, ~] = eigs(L, Z, k, 'SM');

end