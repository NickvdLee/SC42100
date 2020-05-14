clear; close all

load -ascii twitter.mat
W = spconvert(twitter);

W = sparse([W zeros(size(W,1), size(W,1)-size(W,2))]);  % Make square
w = W*ones(size(W,2),1);
D = diag(w);
P = pinv(D)*W;


%% PageRank
T = 1:100;
mu = ones(size(P,1),1);
Beta = 0.15;

y = zeros(size(W,1),size(T,2));
for t = T
    y(:,T+1) = (1 - Beta).*P'*y(:,T)+ Beta.*mu;
end

[PageRank, Indices] = maxk(y(:,end),5);


