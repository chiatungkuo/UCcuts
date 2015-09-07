function u = unifiedcut(La, mincuts, d, alpha, tol, q)
% 
% Find a unified cut from a list of normalized Laplacians
% 
% Inputs
%   La: cell array containing a list of normalized Laplacians
%   mincuts: columnwise-stacked min-cuts for each individual graph in the same order as La
%   d: average node agrees from the original graphs (column vector)
%   alpha: weighting parameter for agreement between output cut to each
%       individual min-cuts
%   tol (optional): a tolerance ratio in (0, 1) for the relaxed orthogonality condition
%   q (optional): number of eigenvectors to look for as candidates
%
% Output
%   u: the unified cut
% 

if nargin < 6, q = 5; end
if nargin < 5, tol = 0.5; end

k = size(La(:), 1);
L = zeros(size(La{1}));
for i = 1:k
    L = L + La{i};
end
Vn = mincuts*mincuts';

A = 1/k*L - alpha*Vn;
A = tril(A) + tril(A, -1)';     % redefine A to ensure it is indeed 'symmetric' to Matlab
[eigvec, ~] = eigs(A, q, 'SA');

for i = 1:q
    u = eigvec(:, i);
    if abs(u'*d.^0.5) < tol*norm(d.^0.5), break; end
end
end
