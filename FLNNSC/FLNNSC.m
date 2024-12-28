function [Z,W,H,obj] = FLNNSC(X,para,maxIter,tol)
% for obtaining similarity matrix
%Input:
%       X:data matrix��size=n*d
%           opts.maxIter: 
%           opts.tol: threshold for judging convergence
%           opts.lambda: 
%           opts.c: 
%Output:
%       Z: representation matrix��size=n*n

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

%     % ��1��У�������ķ���
%     % compute objective function 
%     obj(iter) = norm(X - X*Z, 'fro')^2;
%     % check convergence
%     if iter > 1 && abs(obj(iter) - obj(iter-1)) < tol
%         break;
%     end

    % ��2��У�������ķ���
    obj(iter) = norm(Z-Z_prev,'fro');
    % if iter > 1 && obj(iter) < tol
    %     break;
    % end
    
    % % ��3��У�������ķ���
    % obj(iter) = loss_value(H,Z,L,W,para.alpha,para.beta);

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

