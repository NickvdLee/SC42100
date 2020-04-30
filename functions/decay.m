function [decays] = decay(G,delta)
n = length(G.Nodes.Name);
dist = distances(G);
decays = zeros(n,1);

for i=1:n
    for j=1:n
        if i==j
           continue
        else
            decays(i) = decays(i) + delta^dist(i,j);
        end
    end
end
end

