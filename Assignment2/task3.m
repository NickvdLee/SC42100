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
alpha = 1e-4;      % step size
iterations = 150;


% save xf over all iterations and initialise. 
xfs = cell(4,iterations);
%xfs(4,:) = {zeros(4,1)};

% Set up lambda
nu = cell(1,iterations);
nu(:,:) = {ones(4,1)};

% save inputs and cost for all nodes
inputs = cell(4,1);
cost = cell(4,1);

figure();
hold on


% Solve the problem in one time step:
for i = 1:iterations
    for j = 1:2
        % Making S and T. X_N+1 = Tx_0 + Su_n
        [T,S,W] = mpc_mtrx(A{j},B{j},Tfinal);
        % x = Tx0 + Su
        % xN = A^N x0 + Wu
        
        % Optimising over the control input. 
        cvx_begin quiet
            variable u(2*Tfinal,1)
            variable xf(4,1)
            if j == 1
                minimize((T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u + nu{1,i}'*xf)
            elseif j == 2
                minimize((T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u - nu{1,i}'*xf)
            end
            subject to
                A{j}^Tfinal*x0{j}+W*u == xf;
                u'*u <= umax^2;
        cvx_end
        
        % write to cell
        xfs(j,i) = {xf};
        
        % save for all nodes
        cost{j} = (T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u;
        inputs{j} = u;
        disp(['i: ',num2str(i),' j: ',num2str(j)])
    end
    
    nu{1,i+1} = nu{1,i} + alpha*(xfs{1,i} - xfs{2,i});
       
    plane1 = plot(i,xfs{1,i}(2),'ro');
    plane2 = plot(i,xfs{2,i}(2),'bo');
    % plane3 = plot(i,xfs{3,i}(2),'ko');
    % plane4 = plot(i,xfs{4,i}(2),'go');
    % err = plot(i,abs(sum(xfs{1,i}) - sum(xfs{2,i})),'go');
    drawnow;
end

