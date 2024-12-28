function [Z1,W,H,lap_matrix] = nonlinear(X,W,H,Z1,para)
% FLNNSCC 的非线性部分
%% 设定参数para
alpha = para.alpha;
beta = para.beta;
lambda = para.lambda;
mu = para.mu;
knn = para.knn;
elpson = para.elpson;
%% 更新W和H
[d, n] = size(X);
for i = 1:n
    mapping_x = trigonometric_poly(X(:,i), d)';  % 将d*1的样本映射为5d*1的维度，phi(x)
    W = W - mu*lambda*((H(:,i)-H*Z1(:,i)).*(1.-tanh(W*mapping_x).^2)*mapping_x'+beta*W);  % 更新W
    H(:,i) = tanh(W*mapping_x);  % 计算H，得到映射后的数据；
end
%% 更新Z1
% calculate Laplacian matrix
[pairs,wcost,numpairs]=get_nn_graph(X,knn);
R = zeros(n,numpairs);
for j=1:numpairs
    R(pairs(1,j)+1,j) = wcost(j);
    R(pairs(2,j)+1,j) = -wcost(j);
end
R = R/(knn-1);
lap_matrix = 0.5*(R*R');
% use lyap function to compute Z1
A = H'*H;
B = alpha*(lap_matrix + elpson*eye(size(lap_matrix)));
C = -H'*H;
Z1 = lyap(A,B,C);
end