function value = loss_value(H,Z,L,W,alpha,beta)
%LOSS_VALUE 此处显示有关此函数的摘要
%   此处显示详细说明
value = 1/2*norm(H-H*Z,'fro')^2 + alpha/2*trace(Z*L*Z') + beta/2*norm(W,'fro')^2;

end

