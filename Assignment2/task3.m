clc; clear all; close all;
%{ 
Written by:
    Pim de Bruin
    Nick van der Lee
%}

%% Load data from given file
% state vector is (x,y,xdot,ydot)
run('aircraft.m');
n = size(A1,1);

A = {A1, A2, A3, A4};
B = {B1, B2, B3, B4};
x0 = {x01, x02, x03, x04};

% Systems all controllable -> stabilizable
% Systems are not stable

%% MPC
alpha = 0.5;      % step size
iterations = 20;


% save xf over all iterations and initialise. 
xfs = cell(1,iterations);
xfs(1,:) = {zeros(4,1)};

% Set up lambda
lambda = cell(4,iterations);
lambda(:,:) = {ones(1,4)};

% save inputs and cost for all nodes
inputs = cell(4,1);
cost = cell(4,1);

% Solve the problem in one time step:
for i = 1:iterations
    for j = 1:4
        % Making S and T. X_N+1 = Tx_0 + Su_n
        [T,S,W] = mpc_mtrx(A{j},B{j},Tfinal);
        % x = Tx0 + Su
        % xN = A^N x0 + Wu

        % Optimising over the control input. 
        cvx_begin quiet
            variable u(2*Tfinal,1)
            variable Lambda(1,4)
            minimize((T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u)% + Lambda*((A{j}^Tfinal*x0{j}+W*u) - xfs{1,i}))
            subject to
                A{j}^Tfinal*x0{j}+W*u == xfs{1,i};
                u'*u <= umax^2;
        cvx_end
        
        % save for all nodes
        cost{j} = (T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u;
        inputs{j} = u;
        lambda{j,i} = Lambda;
        disp(['i: ',num2str(i),' j: ',num2str(j)])
    end
    xfs{1,i+1} = xfs{1,i} - alpha*(lambda{1,i}' - lambda{2,i}');      % update xf based on subgradient. 
end

