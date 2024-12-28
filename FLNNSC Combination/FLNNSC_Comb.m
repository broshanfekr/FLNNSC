function [Z,W,H,obj] = FLNNSC_Comb(X,para,maxIter,tol)
% FLNNSC的Combination形式
%Input:
%       X: data matrix，size=n*d
%       para:
%           para.alpha: 
%           para.beta:
%           para.lambda: 
%           para.mu: 
%           para.knn: 
%           para.elpson:
%Output:
%       Z: representation matrix，size=n*n

%% 初始化
[d,n] = size(X);
W = eye(5*d);  % the mapped input = d(1+2*order)
H = zeros(5*d,n);
Z = zeros(n,n);
Z1 = zeros(n,n);
%% 迭代计算表示矩阵Z
[Z2] = linear(X,para.alpha,para.knn,para.elpson);  % 计算Z2
obj = zeros(maxIter, 1);
for iter = 1:maxIter 
    Z_prev = Z;
    [Z1,W,H,~] = nonlinear(X,W,H,Z1,para);  % 更新Z1
    Z = para.lambda*Z1 + (1-para.lambda)*Z2;  % 更新Z

%     % 第1种校验收敛的方法
%     % compute objective function 
%     obj(iter) = norm(X - X*Z, 'fro')^2;
%     % check convergence
%     if iter > 1 && abs(obj(iter) - obj(iter-1)) < tol
%         break;
%     end

    % 第2种校验收敛的方法
    obj(iter) = norm(Z-Z_prev,'fro')^2;
    if iter > 1 && obj(iter) < tol
        break;
    end
    
%     % 第3种校验收敛的方法
%     obj(iter) = loss_value(H,Z,L,W,para.alpha,para.beta);

%     % 第4种校验收敛的方法
%     obj(iter) = norm(Z-Z_prev,'Inf');
%     if iter > 1 && obj(iter) < tol
%         break;
%     end

%     % 第5种校验收敛的方法
%     obj(iter) = max(max(abs(Z-Z_prev)));
%     if iter > 1 && obj(iter) < tol
%         break;
%     end

end
fprintf('迭代次数：%d',iter);
end

