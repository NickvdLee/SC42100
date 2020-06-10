function [T,S,W] = mpc_mtrx(A,B,N)
% A, B = ss matrices
% N = Np = Nc time horizons
% t = no. of iterations until Tf

C = eye(size(A)); % Full state feedback

nx = size(A,2); nu = size(B,2); ny = size(C,1);

%% Exogenous response, T: xa = Tx0
T = zeros(ny*N,nx);
for k=0:N-1
    T(k*ny+1:(k+1)*ny,:) = A^k;
end

%% Input response, S: xi = Su
S = zeros(ny*N,nu*N);
for k=1:N-1 % To N-1
    for i=0:k-1
        S(k*ny+1:(k+1)*ny,i*nu+1:(i+1)*nu) = A^(k-1-i)*B;
    end
end
% State sequence prediction: x = Tx0 + Su
%% Final state prediction W: A^N x0 + Wu = x(N)
W = zeros(nx,nu*N);
for i=0:N-1 % To N
    W(1:nx,i*nu+1:(i+1)*nu) = A^(N-i-1)*B;
end

end