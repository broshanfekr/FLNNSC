function [Z,W,H,lap_matrix] = flnnsc_paraUpdate(X, W, H, Z, para)
%% parameter update
% parameters include \mu,\alpha,\beta, knn, \elpson
mu = para.mu;
alpha = para.alpha;
beta = para.beta;
knn = para.knn;
elpson = para.elpson;
%% 更新W和H
[d, n] = size(X);
for i = 1:n
    mapping_x = trigonometric_poly(X(:,i), d)';  % 将d*1的样本映射为5d*1的维度，phi(x)
    W = W - mu*((H(:,i)-H*Z(:,i)).*(1.-tanh(W*mapping_x).^2)*mapping_x'+beta*W);  % 更新W
    H(:,i) = tanh(W*mapping_x);  % 计算H，得到映射后的数据；
end
%% 更新Z
% calculate Laplacian matrix
[pairs,wcost,numpairs]=get_nn_graph(X,knn);
R = zeros(n,numpairs);
for j=1:numpairs
    R(pairs(1,j)+1,j) = wcost(j);
    R(pairs(2,j)+1,j) = -wcost(j);
end
R = R/(knn-1);
lap_matrix = 0.5*(R*R');

% % construct similarity matrix S
% S = zeros(n, n);
% sigma = 1;
% for i = 1:n
%     for j = i+1:n
%         S(i, j) = exp(-norm(X(:, i) - X(:, j))^2 / (2*sigma^2));
%         S(j, i) = S(i, j);
%     end
% end
% % compute degree matrix D and Laplacian matrix
% D = diag(sum(S, 2));
% lap_matrix = D - S;


% use lyap function to compute Z
A = H'*H;
B = alpha*(lap_matrix + elpson*eye(size(lap_matrix)));
C = -H'*H;
Z = lyap(A,B,C);

end