function [Z,W,H,obj] = FLNNSC(X,para,maxIter,tol)
% for obtaining similarity matrix
%Input:
%       X:data matrix，size=n*d
%           opts.maxIter: 
%           opts.tol: threshold for judging convergence
%           opts.lambda: 
%           opts.c: 
%Output:
%       Z: representation matrix，size=n*n

%% initialization
[d,n] = size(X);
W = eye(5*d); % the mapped input = d(1+2*order)
H = zeros(5*d,n);
Z = zeros(n,n);

%% Iterative calculation of the representation matrix Z
obj = zeros(maxIter, 1);  % 1.3 Verify convergence needs
for iter = 1:maxIter 
    % update Z
    Z_prev = Z;
    [Z,W,H,L] = flnnsc_paraUpdate(X, W, H, Z, para);

%     % 第1种校验收敛的方法
%     % compute objective function 
%     obj(iter) = norm(X - X*Z, 'fro')^2;
%     % check convergence
%     if iter > 1 && abs(obj(iter) - obj(iter-1)) < tol
%         break;
%     end

    % 第2种校验收敛的方法
    obj(iter) = norm(Z-Z_prev,'fro');
    % if iter > 1 && obj(iter) < tol
    %     break;
    % end
    
    % % 第3种校验收敛的方法
    % obj(iter) = loss_value(H,Z,L,W,para.alpha,para.beta);

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

