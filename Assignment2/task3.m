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
iterations = 500;

% save xf over all iterations and initialise. 
xfs = cell(4,iterations);
xfs(:,:) = {zeros(4,1)};

%%
% Set up lambda
lambda = cell(n,n,iterations);
lambda(:,:,:) = {ones(4,1)};
for k=1:n
    lambda(k,k,:) = {zeros(4,1)};
end

% save inputs and cost for all nodes
inputs = cell(n,1);
cost = cell(n,1);
%%
% Solve the problem in one time step:
T = cell(1,n);
S = cell(1,n);
W = cell(1,n);
for i = 1:n
    [t,s,w] = mpc_mtrx(A{i},B{i},Tfinal);
    T{i} = t;
    S{i} = s;
    W{i} = w;
end
for k = 1:iterations
    for i = 1:n
        % Making S and T. X_N+1 = Tx_0 + Su_n
%         [T,S,W] = mpc_mtrx(A{i},B{i},Tfinal);
        % x = Tx0 + Su
        % xN = A^N x0 + Wu
        
        % Optimising over the control input. 
        cvx_begin quiet
            variable u(2*Tfinal,1)
            x = T{i}*x0{i} + S{i}*u;
            xf = A{i}^Tfinal*x0{i}+W{i}*u;
            lij = 0;            
            lji = 0;
            for j=1:n
                lij = lij + lambda{i,j,k};
                lji = lji + lambda{j,i,k};
            end
            theta = repmat(xf,[Tfinal,1]);
            minimize((x-theta)'*(x-theta) + u'*u + lij'*xf - lji'*xf)
            subject to
                u'*u <= umax^2;
        cvx_end
        
        inputs{i} = u;
        xfs{i,k} = xf;
        disp(['i: ',num2str(k),' j: ',num2str(i)])
    end
    
    % Update lij and lji
    for i=1:n
        for j=1:n
            if i==j
                continue
            end
            lambda{i,j,k+1} = lambda{i,j,k} + alpha*(xfs{i,k}-xfs{j,k});
            lambda{j,i,k+1} = lambda{j,i,k} + alpha*(xfs{j,k}-xfs{i,k});
        end
    end

    disp('Final states:');
    disp([xfs{1,k} xfs{2,k} xfs{3,k} xfs{4,k}]);
    disp(norm(diff([xfs{1,k} xfs{2,k} xfs{3,k} xfs{4,k}],n-1,2)));
end

