function [c] = closeness(G)
n = length(G.Nodes.Name);
c = n./sum(distances(G),2);
end

