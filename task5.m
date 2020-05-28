%{
task 5.

Written by:
Pim de Bruin        4545702
Nick van der Lee    4144600
%}


close all; clear all; clc;
addpath('functions/');
% addpath('cvx/');
% cvx_startup
load -ascii traffic.mat 
load -ascii traveltime.mat
load -ascii capacities.mat
load -ascii flow.mat


%% Setup graph G
% get edges from B
Links = zeros(size(traffic,2),2);
for i = 1:size(traffic,2)
    Links(i,1) = find(traffic(:,i) == 1);
    Links(i,2) = find(traffic(:,i) == -1);
end

% also minimum travel time
% Links(:,3) = traveltime;

% Make graph using empty network weights
G = digraph(Links(:,1)', Links(:,2)', traveltime');

%% 5.a
[p1_13, d1_13, e1_13]  = shortestpath(G,1,13);

%% 5.b
% matlab will treat the weight as the capacity of the edge. Therefore the
% weight vector needs to be swapped out for the capacity vector. 

Links(:,3) = capacities;
G = digraph(Links(:,1)',Links(:,2), Links(:,3)');
[mf1_13,~] = maxflow(G,1,13);

%% 5.c
nu = traffic*flow; % nu = B*f

%% 5.d
% Define important parameters. 
B = traffic;
l = traveltime;
M = size(traffic,2);
c = capacities;
lambda = [nu(1) zeros(1,size(traffic,1)-1)]';
mu = [zeros(1,size(traffic,1)-1) nu(1)]';
One = ones(M,1);


% Compute social optimum
cvx_begin quiet
    variable f(M)
    minimize sum((l.*c).*inv_pos(One-f.*(1.*inv_pos(c))))   % Use correct cost function
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fstar = f; % The social optimum

% Plot the flow over the graph
figure()
Gplot = digraph(Links(:,1)', Links(:,2)', fstar');
W = Gplot.Edges.Weight./max(Gplot.Edges.Weight);
plot(Gplot,'k','EdgeLabel',round(Gplot.Edges.Weight),'LineWidth',3*W,'ArrowSize',15*W,'Layout','force');

%% 5.e
% Compute wardrup optimum using the appropriate cost function
cvx_begin quiet
    variable f(M)
    minimize sum(-l.*c.*log(c-f))
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fw = f;

% Plot the flow over the graph. 
figure()
Gplot = digraph(Links(:,1)', Links(:,2)', fw');
W = Gplot.Edges.Weight./max(Gplot.Edges.Weight);
plot(Gplot,'k','EdgeLabel',round(Gplot.Edges.Weight),'LineWidth',3*W,'ArrowSize',15*W,'Layout','force');
fo = f; % Wardrup optimum


%% 5.f
% Compute toll vector w
dacc = c.*l./((c-fstar).*(c-fstar));
w = fstar.*dacc;

cvx_begin quiet
    variable f(M)
    minimize sum(-l.*c.*log(c-f)+w.*f)
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fw = f; % wardrup optimal vector using tolls


%% 5.g
% compute the new fstar using the new cost function:
cvx_begin
    variable f(M)
    minimize sum((l.*c).*inv_pos(One-f.*(1.*inv_pos(c)))-f.*l)
    subject to
        B*f == lambda - mu
        0 <= f <= c
cvx_end
fstar2 = f; % new fstar using new cost function


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
fw2 = f; % new wardrup optimum using new cost function and new tolls. 


% Plot added travel times
figure()
Gplot = digraph(Links(:,1)', Links(:,2)', (l./(1-fstar2./c) - traveltime)');
W = Gplot.Edges.Weight./max(Gplot.Edges.Weight);
h = plot(Gplot,'k','EdgeLabel',round(Gplot.Edges.Weight,2),'ArrowSize',0,'LineWidth',3*W,'Layout','force');