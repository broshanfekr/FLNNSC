function [Z2] = linear(X,alpha,knn,elpson)
% FLNNSCC的线性部分
%% 计算Z2
% calculate Laplacian matrix
[~,n] = size(X);
[pairs,wcost,numpairs]=get_nn_graph(X,knn);
R = zeros(n,numpairs);
for j=1:numpairs
    R(pairs(1,j)+1,j) = wcost(j);
    R(pairs(2,j)+1,j) = -wcost(j);
end
R = R/(knn-1);
lap_matrix = 0.5*(R*R');
% use lyap function to compute Z2
A = X'*X;
B = alpha*(lap_matrix + elpson*eye(size(lap_matrix)));
C = -X'*X;
Z2 = lyap(A,B,C);
end

