function [c] = closeness(G)
% IN: GRAPH
% OUT: SCALAR R^numnodes(G)

n = length(G.Nodes.Name);
c = n./sum(distances(G),2);
end

