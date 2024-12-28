function [Z,W,H,obj] = FLNNSC_Comb(X,para,maxIter,tol)
% FLNNSC��Combination��ʽ
%Input:
%       X: data matrix��size=n*d
%       para:
%           para.alpha: 
%           para.beta:
%           para.lambda: 
%           para.mu: 
%           para.knn: 
%           para.elpson:
%Output:
%       Z: representation matrix��size=n*n

%% ��ʼ��
[d,n] = size(X);
W = eye(5*d);  % the mapped input = d(1+2*order)
H = zeros(5*d,n);
Z = zeros(n,n);
Z1 = zeros(n,n);
%% ���������ʾ����Z
[Z2] = linear(X,para.alpha,para.knn,para.elpson);  % ����Z2
obj = zeros(maxIter, 1);
for iter = 1:maxIter 
    Z_prev = Z;
    [Z1,W,H,~] = nonlinear(X,W,H,Z1,para);  % ����Z1
    Z = para.lambda*Z1 + (1-para.lambda)*Z2;  % ����Z

%     % ��1��У�������ķ���
%     % compute objective function 
%     obj(iter) = norm(X - X*Z, 'fro')^2;
%     % check convergence
%     if iter > 1 && abs(obj(iter) - obj(iter-1)) < tol
%         break;
%     end

    % ��2��У�������ķ���
    obj(iter) = norm(Z-Z_prev,'fro')^2;
    if iter > 1 && obj(iter) < tol
        break;
    end
    
%     % ��3��У�������ķ���
%     obj(iter) = loss_value(H,Z,L,W,para.alpha,para.beta);

%     % ��4��У�������ķ���
%     obj(iter) = norm(Z-Z_prev,'Inf');
%     if iter > 1 && obj(iter) < tol
%         break;
%     end

%     % ��5��У�������ķ���
%     obj(iter) = max(max(abs(Z-Z_prev)));
%     if iter > 1 && obj(iter) < tol
%         break;
%     end

end
fprintf('����������%d',iter);
end

