close all; clear all; clc;
addpath('functions/');

%% Construct graph
cities = {'Alkmaar','Haarlem','Leiden','Den Haag','Hoorn','Amsterdam','Schiphol','Rotterdam','Dordrecht','Utrecht','Almere','Lelystad','Amersfoort'};
n = length(cities);

W = zeros(n);
connections = {[2 5 6],...      % Alkmaar
               [1 6 3],...      % Haarlem
               [2 7 4],...      % Leiden
               [3 8],...        % Den Haag
               [1 6],...        % Hoorn
               [1 2 5 7 11],... % Amsterdam
               [3 6 8 10],...   % Schiphol
               [4 7 9 10],...   % Rotterdam
                8,...           % Dordrecht
               [7 3 8 13],...   % Utrecht
               [6 12],...       % Almere
                11,...          % Lelystad
               [6 10]};         % Amersfoort
for i = 1:length(connections)
    idxs = connections{i};
    for j = 1:length(idxs)
        W = connect(W,i,idxs(j));
    end
end
G = graph(W,cities);

%% Centrality measures
c_bona = bonacich(G);
c_close = closeness(G);

delta = [0.25 0.5 0.75];
c_decay = zeros(n,length(delta));
for d=1:length(delta)
    c_decay(:,d) = decay(G,delta(d));
end

%% c_betw = centrality(G,'betweenness');
c_betw = zeros(n,1);
V = G.Nodes;
for s=1:numnodes(G) % Represent node as numeric
    S = []; % Empty stack
    P = cell(numnodes(G),1); % Empty list?
    sigma = zeros(numnodes(G),1);
    sigma(s) = 1;
    d = -ones(numnodes(G),1);
    d(s) = 0;
    Q = [s]; % Queue
    while ~isempty(Q)
        v = Q(1); Q = Q(2:end);
        S = [S v];
        
        N = neighbors(G,v);
        for k=1:length(N)
            w = N(k);
            % w found for the first time?
            if d(w) < 0
                Q = [Q w];
                d(w) = d(v) + 1;
            end
            % shortest path to w via v?
            if d(w) == d(v) + 1
                sigma(w) = sigma(w) + sigma(v);
                if isempty(P{w})
                    P{w} = [v];
                else
                    P{w} = [P{w} v];
                end
            end
        end
    end
    
    delta = zeros(numnodes(G),1);
    while ~isempty(S)
        w = S(1); S = S(2:end);
        for k=1:length(P{w})
            delta(v) = delta(v) + sigma(v)/sigma(w)*(1+delta(w));
        end
        if w ~= s
            c_betw(w) = c_betw(w)+delta(w);
        end
    end
    
end

%% Prints
enum = [findnode(G,'Amsterdam'),...
        findnode(G,'Schiphol'),...
        findnode(G,'Rotterdam'),...
        findnode(G,'Utrecht')];
dataTable = table(cities(enum)',c_bona(enum),c_close(enum),c_decay(enum,1),c_decay(enum,2),c_decay(enum,3),c_betw(enum), ...
    'VariableNames',{'City','Bonacich','Closeness','Decay (0.25)','Decay (0.5)','Decay (0.75)','Betweenness'});
disp(dataTable)

function W = connect(W,a,b)
    W(a,b) = 1;
    W(b,a) = 1;
end