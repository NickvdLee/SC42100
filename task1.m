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
c_decay = zeros(n,0);
for d=delta
    c_decay(:,end+1) = decay(G,d);
end

c_betw = betweenness(G);

%% Prints
enum = [findnode(G,'Amsterdam'),...
        findnode(G,'Schiphol'),...
        findnode(G,'Rotterdam'),...
        findnode(G,'Utrecht')];
dataTable = table(cities(enum)',c_bona(enum),c_close(enum),c_decay(enum,1),c_decay(enum,2),c_decay(enum,3),c_betw(enum), ...
    'VariableNames',{'City','Bonacich','Closeness','Decay (0.25)','Decay (0.5)','Decay (0.75)','Betweenness'});
disp(dataTable)

%% Graph construction functions
function W = connect(W,a,b)
    W(a,b) = 1;
    W(b,a) = 1;
end