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
% Generate first xf for case of zero input
xf = zeros(4,1);

inputs = cell(4,1);
cost = cell(4,1);

% Solve the problem in one time step:

for j = 1:4
    % Making S and T. X_N+1 = Tx_0 + Su_n
    [T,S,W] = mpc_mtrx(A{j},B{j},Tfinal);
    % x = Tx0 + Su
    % xN = A^N x0 + Wu
    
    % Optimising over the control input. 
    cvx_begin quiet
    variable u(2*Tfinal,1)
    minimize((T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u)
    subject to
        A{j}^Tfinal*x0{j}+W*u == xf;
        u'*u <= umax^2;
    cvx_end
    
    
    cost{j} = (T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u;
    inputs{j} = u;
end
