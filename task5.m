close all; clear all; clc;
addpath('functions/');
addpath('cvx/');
cvx_startup
load -ascii traffic.mat 
load -ascii traveltime.mat
load -ascii capacities.mat
load -ascii flow.mat


%% Setup graph G
% get edges from B
Links = zeros(size(traffic,2),3);
for i = 1:size(traffic,2)
    Links(i,1) = find(traffic(:,i) == 1);
    Links(i,2) = find(traffic(:,i) == -1);
end

% also minimum travel time
Links(:,3) = traveltime;

% Make graph using empty network weights
G = digraph(Links(:,1)', Links(:,2)', Links(:,3)');

%% 5.a
[p1_13, d1_13, e1_13]  = shortestpath(G,1,13);

%% 5.b
% matlab will treat the weight as the capacity of the edge. Therefore the
% weight vector needs to be swapped out for the capacity vector. 

Links(:,3) = capacities;
G = digraph(Links(:,1)',Links(:,2), Links(:,3)');
[mf1_13,~] = maxflow(G,1,13);

%% 5.c
nu = traffic*flow;

%% 5.d
B = traffic;
l = traveltime;
M = size(traffic,2);
c = capacities;
lambda = [nu(1) zeros(1,size(traffic,1)-1)]';
mu = [zeros(1,size(traffic,1)-1) nu(1)]';
One = ones(M,1);

cvx_begin
    variable f(M)
    minimize sum((l.*c).*inv_pos(One-f.*(1.*inv_pos(c))))
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fstar = f;


%% 5.e
cvx_begin
    variable f(M)
    minimize sum(-l.*c.*log(c-f)/log(exp(1)))
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fo = f;

%% 5.f
dacc = c.*l./((c-fstar).*(c-fstar));
w = fstar.*dacc;

cvx_begin
    variable f(M)
    minimize sum(-l.*c.*log(c-f)+w.*f)
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fw = f;

%% 5.g
% compute the new fstar using the new cost function:
cvx_begin
    variable f(M)
    minimize sum((l.*c).*inv_pos(One-f.*(1.*inv_pos(c)))-f.*l)
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fstar2 = f;

% compute the wardrup equilibrium using the tolls:
dacc = c.*l./((c-fstar2).*(c-fstar2));
w2 = fstar2.*dacc - l; % -l due to the altered cost function, so c'is different. 

cvx_begin
    variable f(M)
    minimize sum(-l.*c.*log(c-f)+w2.*f)
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fw2 = f;
