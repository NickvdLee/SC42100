clc; clear all; close all;
set(0,'defaultTextInterpreter','latex');
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
alpha = .15;      % step size
iter = 100;

% save xf over all iterations and initialise. 
xfs = cell(4,iter);
xfs(:,:) = {zeros(4,1)};
xs = cell(4,iter);
xs(:,:) = {zeros((Tfinal+1)*size(A1,1),1)};

%%
% Set up lambda
lambda = cell(n,n,iter);
lambda(:,:,:) = {ones(4,1)};
for r=1:n
    lambda(r,r,:) = {zeros(4,1)};
end

% save inputs and cost for all nodes
inputs = cell(n,iter);
cost = cell(n,1);

% Solve the problem in one time step:
T = cell(1,n); S = cell(1,n); W = cell(1,n);
for i = 1:n
    [t,s,w] = mpc_mtrx(A{i},B{i},Tfinal);
    T{i} = sparse(t); S{i} = sparse(s); W{i} = w;
end

for r = 1:iter
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
                lij = lij + lambda{i,j,r};
                lji = lji + lambda{j,i,r};
            end
            theta = repmat(xf,[Tfinal,1]);
            minimize(x'*x + xf'*xf + u'*u + lij'*xf - lji'*xf)
            subject to
                u'*u <= umax^2;
                xf >= zeros(4,1);
        cvx_end
        
        inputs{i,r} = u;
        xfs{i,r} = xf;
        xs{i,r} = [T{i}*x0{i}+S{i}*u;xf];
        disp(['i: ',num2str(r),' j: ',num2str(i)])
    end
    
    % Update lij and lji
    for i=1:n
        for j=1:n
            if i==j
                continue
            end
            lambda{i,j,r+1} = lambda{i,j,r} + alpha*(xfs{i,r}-xfs{j,r});
            lambda{j,i,r+1} = lambda{j,i,r} + alpha*(xfs{j,r}-xfs{i,r});
        end
    end

    disp('Final states:');
    disp([xfs{1,r} xfs{2,r} xfs{3,r} xfs{4,r}]);
    disp('Convergence:');
    disp(norm(diff([xfs{1,r} xfs{2,r} xfs{3,r} xfs{4,r}],n-1,2)));
end

%% Plots

% Xs
figure(1)
colors = {[0 0 1],[1 0 0],[0 1 0],[0 1 1]};
for r = 1:iter
    for i=1:4
        X = reshape(xs{i,r},[4 6]);
        plot(X(1,:),X(2,:),'-','Color',(colors{i}*r + [.75 .75 .75].*(iter-r))/iter);
        hold on
        axis(([X(1,6)-1 X(1,6)+1 X(2,6)-1 X(2,6)+1]*r + [-10 10 -10 10]*(iter-r))/iter);
    end
    hold off
    drawnow
    pause(0.1);  
end


% Xfs
Xf = zeros(4,iter,4);
p0 = gobjects(4);
shapes = {'-','--',':','-.'};
markers = {'o','x','^','s'};

for r = 1:iter
    for j = 1:n     % 4 planes
        for i =1:4  % 4 states
            Xf(j,r,i) = xfs{j,r}(i);
        end
    end
end

figure(2)
hold on;
for i = 1:4         % 4 states
    for j = 1:n     % 4 planes
        p0(i,j) = plot(1:iter,Xf(j,:,i),shapes{j},'Color',colors{i});
    end
end
xlabel('Iterations')
ylabel('State value')
ylim([0 4]);
title("Private $$x_f$$ for all 4 planes over all iterations")

p1 = plot(1:iter,Xf(1,:,1),shapes{1},'Color',colors{1}); % plot first one again for the legend. (otherwise it has 2 legend entries)
lgd = legend([p1 p0(1,2:end),p0(:,1)'],'Plane 1','Plane 2','Plane 3','Plane 4','x','y','xdot','ydot');
lgd.NumColumns = 2;


% inputs
U = zeros(n,iter);
for r = 1:iter
    for j = 1:n     % 4 planes
        U(j,r) = inputs{j,r}'*inputs{j,r};
    end
end

figure(3);
hold on;
p2 = zeros(4,1);
for j = 1:n
    p2(j) =  plot(1:iter,U(j,:),shapes{j},'Color','k');   
end
legend('Plane 1','Plane 2','Plane 3','Plane 4','Location','northeast');
title('$$u^Tu$$ per plane per time step');
xlabel('Iterations')
ylabel('$$u^Tu$$')


% states
X = zeros(4,Tfinal+1,n); % states X time steps X planes
U = zeros(2,Tfinal,n);   % inputs X time steps X planes 
p4 = gobjects(4,1);
p5 = gobjects(4,1);


figure(4);
sgtitle('States of all 4 planes');
subplot(3,2,[1,2,3,4])
    hold on;
    for j = 1:n     % 4 planes
        X(:,:,j) = reshape(xs{j,iter},[4,6]);
        plot(X(1,1:end-1,j),X(2,1:end-1,j),strcat(shapes{j},markers{j}),'Color','k');
        hold on;
    end
    xlabel('x position')
    ylabel('y position')
    legend('Plane 1','Plane 2','Plane 3','Plane 4','Location','east');
subplot(3,2,[5,6])
    hold on;
    for j = 1:n     % 4 planes
        p4(j) = plot(1:Tfinal,X(3,1:end-1,j),shapes{j},'Color',colors{3});
        p5(j) = plot(1:Tfinal,X(4,1:end-1,j),shapes{j},'Color',colors{4});
    end
    lgd = legend([p4(1) p5(1)],'xdot','ydot');
    xlabel('Time steps')
    ylabel('velocity')








