function u = consensus(cuts)
%
% Compute a consensus cut from a given set of cuts.
% The method is adapted from "Consensus clustering in complex networks" by 
% Lancichinetti and Fortunato in 2012.
%
% Input
%   cuts: n by m matrix of m cuts of size n
%
% Output
%   u: the consensus cut
%

[n, m] = size(cuts);
% consensus matrix recording the proportion of agreement for each node pair
C = zeros(n, n);    

for i = 1:n
    for j = i:n
        for k = 1:m
            cut = cuts(:, k);
            if cut(i)*cut(j) > 0, C(i, j) = C(i, j) + 1; end
        end
    end
end
A = (triu(C, 1)' + C)/n;
D = diag(sum(A, 2));
L = eye(n) - D^(-1/2)*A*D^(-1/2);
[V, ~] = eigs(L, 3, 'SM');
u = V(:, 2);

end