function u = contrastcut(La, Lb, d, beta, tol, q)
% 
% Find a contrast cut between two sets of graphs (Laplacians)
% 
% Inputs
%   La: cell array containing normalized Laplacians for which the contrast 
%       with low cost is desired
%   Lb: cell array containing normalized Laplacians for which the contrast
%       cut gives high cost
%   d: average node agrees from the original graphs of La (column vector)
%   beta: weighting parameter for Lb
%   tol (optional): a tolerance ratio in (0, 1) for the relaxed orthogonality condition
%   q (optional): number of eigenvectors to look for as candidates
%
% Output
%   u: the contrast cut
%

if nargin < 6, q = 5; end
if nargin < 5, tol = 0.5; end

ka = size(La(:), 1);
kb = size(Lb(:), 1);

A = zeros(size(La{1}));
for i = 1:ka
    A = A + La{i};
end

B = zeros(size(Lb{1}));
for j = 1:kb
    B = B + Lb{j};
end

Z = 1/ka*A - beta*1/kb*B;
Z = tril(Z) + tril(Z, -1)'; % make sure it is "numerically" symmetric
[eigvec ~] = eigs(Z, q, 'SA');

for i = 1:q
    u = eigvec(:, i);
    if (abs(u'*d.^0.5) <= tol*norm(d.^0.5)), break; end
end

end
