close all; clear all; clc;
addpath('functions/');
load -ascii twitter.mat
W = spconvert(twitter); clear twitter;

%% Construct graph
l = length(W);
W(l,l) = 0;

% l = 100;
W = W(1:l,1:l);

G = digraph(W);

P = diag(1./sum(W,2))*W; % Fast inverse

mu = ones(l,1);
beta = 0.15; % Typical value

PageRank = mu; % Typical initialization?
M = (1-beta)*P'; % Cache these
N = beta*mu;
for i=1:10
    PageRank = M*PageRank + N;
end
% Pagerank = (eye(l)-M)\N; SS solution

[PR, idxs] = maxk(PageRank,5);

%% Stubborn nodes
% Let's find out if there already are stubborn nodes (sinks)
sources = [];
sinks = [];
in = indegree(G);
out = outdegree(G);
for k=1:numnodes(G)
    if in(k) == 0
        sources = [sources;k];
    elseif out(k) == 0
        sinks = [sinks;k];
    end
end