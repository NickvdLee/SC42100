%{
task 3.

Written by:
Pim de Bruin        4545702
Nick van der Lee    4144600
%}

close all; clear all; clc;
addpath('functions/');
load -ascii twitter.mat
W = spconvert(twitter); clear twitter;

%% Construct graph
l = length(W);
W(l,l) = 0; % Resize by adding proper collums

G = digraph(W);
P = diag(1./sum(W,2))*W; % Fast inverse

mu = ones(l,1);
beta = 0.15; % Typical value


% set up pagerank iteratively
PageRank = mu; 
M = (1-beta)*P';
N = beta*mu;
for i=1:75
    PageRank = M*PageRank + N;
end

[PR, idxs] = maxk(PageRank,5); % find maximum 5 page ranks. 

%% Stubborn nodes
% Let's find out if there already are stubborn nodes (sinks)
in = find(indegree(G) == 0);
out = find(outdegree(G) == 0);

% 1 is already a sink, and the pagerank idxs were
% 2 1 112 9 26 so lets also chose one of them
c = 112;
S = [1 c];
s = 2;
R = [2:c-1 c+1:l];
r = l-s;

% Remove outgoing edges
for edge=outedges(G,c)
    G = rmedge(G,edge);
end

% re compute graph
W = adjacency(G);
P = diag(1./sum(W,2))*W; % Fast inverse
Q = P(R,R);
B = P(R,S);

u = [1;0]; 

% Dynamics
D = distances(G);
D(isinf(D)) = -1;
D = max(max(D));

% Simulate system
T = 150;
y = zeros(r,T);
y(:,1) = rand(r,1);
for t=1:T-1
    y(:,t+1) = Q*y(:,t) + B*u;
end

% Plot
num_plots = 5;
yidxs = sort(randi(r,1,num_plots));
yidxnames = cell(1,num_plots);
for k=1:num_plots
    if yidxs(k) < c
        yidxnames{k} = sprintf('%d',yidxs(k)+1);
    else
        yidxnames{k} = sprintf('%d',yidxs(k)+2);
    end
end

figure(1)
hold off
for i=yidxs
    stairs(1:T,y(i,:))
    hold on
end
legend(yidxnames)
xlabel('k')
ylabel('Opinion')

%% Plot opinion distribution
% Opinions are a continuous number in (0,1)
% Therefore, boxplot!
Y = [];
% Node strength:
% strongstrong strongweak weakstrong weakweak
C = [[1;112] [2;randi(l,1,1)] [randi(l,1,1);26] sort(randi(l,2,1))];
for S=C
    G = digraph(W);
    s = 2;
    R = 1:l;
    R(S) = [];
    r = l-s;
    % Remove outgoing edges
    for i = 1:2
        for edge=outedges(G,S(i))
            G = rmedge(G,edge);
        end
    end
    W = adjacency(G);
    P = diag(1./sum(W,2))*W; % Fast inverse
    Q = P(R,R);
    B = P(R,S);

    u = [1;0];

    % Dynamics
    T = 150;
    y = rand(r,1);
    y(1) = 0.5;
    
    for t=1:T-1
        y = Q*y + B*u;
    end
    Y = [Y y];
end


% Plot bpxplot
figure(2)
boxplot(Y,'Colors','k')

names = cell(1,4);
for k=1:4
    names{k} = sprintf('%d = 1, %d = 0',C(:,k));
end
xticklabels(names)