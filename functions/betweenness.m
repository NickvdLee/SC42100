function bc = betweenness( G )

bc = zeros(numnodes(G),1);
n = numnodes(G);
for s=1:n
    S = [];         % Empty stack
    P = cell(n,1);  % Empty list
    sigma = zeros(n,1); sigma(s) = 1;
    d = -ones(n,1); d(s) = 0;
    Q = Queue;      % Empty queue
    Q.enqueue(s);
    while ~isempty(Q)
        v = Q.dequeue;
        S = [S v];
        N = neighbors(G,v);
        for k=1:length(N)
            w = N(k);
            if d(w)<0 % First visit
                Q.enqueue(w);
                d(w) = d(v)+1;
            end
            if d(w) == d(v)+1 % Shortest path to w via v
                sigma(w) = sigma(w)+sigma(v);
                P{w,1} = [P{w,1} v];
            end                
        end
    end
    delta = zeros(n,1);
    while ~isempty(S)
        w = S(end); S = S(1:end-1);
        for k=1:length(P{w,1})
            v = P{w,1}(k);
            delta(v) = delta(v)+sigma(v)/sigma(w)*(1+delta(w));
        end
        if w~=s
            bc(w) = bc(w) + delta(w);
        end
    end
end
bc = bc./(n^2);