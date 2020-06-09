clc; clear all; close all;
%{ 
Written by:
    Pim de Bruin
    Nick van der Lee
%}

% state vector is (x,y,xdot,ydot)

A1=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B1=[2 0;0 2;3 0;0 3];
x01=[-10;10;-1;1];

A2=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B2=[3 0; 0 3; 7 0; 0 7];
x02=[10;10;1;1];

A3=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B3=[1 0; 0 1; 1.1 0; 0 1.1];
x03=[10;-10;1;-1];

A4=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B4=[6 0;0 6;20 0; 0 20];
x04=[-10;-10;-1;-1];

Tfinal=5;   % randez-vouz time
umax=100;   % Limit on controls
n = 4;      % # of states


% put in cells for loopint
A = {A1, A2, A3, A4};
B = {B1, B2, B3, B4};
x0 = {x01, x02, x03, x04};





%% in the master optimisation loop!
xf = [1; 1; 1; 1];
inputs = cell(4,1);
cost = cell(4,1);

% do for every plane:
for j = 1:4
    % Making S and T. X_N+1 = Tx_0 + Su_n
    T = zeros(n*(Tfinal),n);
    T(1:n,1:n) = eye(n);
    S = zeros(n*(Tfinal),2*(Tfinal));
    for i = 2:Tfinal
        T((n*i)-n+1:n*i,:) = A{j}^(i-1);
        for k = 1:i-1
            S((n*i)-n+1:n*i,2*k-1:2*k) = A{j}^(i-k-1)*B{j};
        end
    end
    
    
    % Optimising over the control input. 
    cvx_begin quiet
    variable u(2*Tfinal,1)
    minimize((T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u)
    subject to
        T(Tfinal*n-n+1:end,:)*x01+S(Tfinal*n-n+1:end,:)*u == xf;
    cvx_end
    
    
    cost{j} = (T*x0{j}+S*u)'*(T*x0{j}+S*u) + u'*u;
    inputs{j} = u;
end
