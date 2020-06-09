function [greek_pi] = bonacich(G)
% IN: GRAPH
% OUT: SCALAR R^numnodes(G)

W = full(adjacency(G));
w = W*ones(size(W,2),1);
P = diag(w)\W;
% Assume there always exists an eigenvalue of P' of 1
greek_pi = null(P'-eye(size(P')));
end
