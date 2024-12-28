function [Aff] = get_Aff(Z,X,aff_type,gamma)
%GET_AFF 从表示矩阵Z得到亲和矩阵A的三种方式
%Input：
%       Z：表示矩阵；
%       X：数据矩阵；
%       aff_type：获得亲和矩阵的方式，有'J1'、'J2'、'J2_nonorm'、'thr'；
%       gamma：其中的参数；如果'J1'，gamma无作用；'J2'、'J2_nonorm'则一般取值[1,10]；
%       'thr'，则取值[0,1]；
%Output：
%       Aff：亲和矩阵；

nX = sqrt(sum(X.^2));

if strcmp(aff_type,'J1')
    Aff =(abs(Z)+abs(Z'))/2;
elseif strcmp(aff_type,'J2')
    Aff = abs(Z'*Z./(nX'*nX)).^gamma;
elseif strcmp(aff_type,'J2_nonorm')
    Aff = abs(Z'*Z).^gamma;
elseif strcmp(aff_type,'thr')
    % 要进行一个阈值操作，然后根据最强的k个关系进行构图，但后续还是'J1'的形式
    Z = thrC(Z,gamma);
    Aff = BuildAdjacency(Z);
end

end


function Cp = thrC(C,ro)
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------
% 一种阈值操作，将输入的表示矩阵C的每列的绝对值最小的几个数置为0，参数ro控制
% 过滤数量，取值范围为[0,1]，为1则不变，越小则过滤得越多，一般取0.4-0.7

if (nargin < 2)
    ro = 1;
end

if (ro < 1)
    N = size(C,2);
    Cp = zeros(N,N);
    [S,Ind] = sort(abs(C),1,'descend');
    for i = 1:N
        cL1 = sum(S(:,i));
        stop = false;
        cSum = 0; t = 0;
        while (~stop)
            t = t + 1;
            cSum = cSum + S(t,i);
            if ( cSum >= ro*cL1 )
                stop = true;
                Cp(Ind(1:t,i),i) = C(Ind(1:t,i),i);
            end
        end
    end
else
    Cp = C;
end
end

function [CKSym,CAbs] = BuildAdjacency(CMat,K)
%--------------------------------------------------------------------------
% This function takes a NxN coefficient matrix and returns a NxN adjacency
% matrix by choosing the K strongest connections in the similarity graph
% CMat: NxN coefficient matrix
% K: number of strongest edges to keep; if K=0 use all the exiting edges
% CKSym: NxN symmetric adjacency matrix
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------
% 源代码存在问题，只是选取了最强的K个关系进行归一化，并没有去弱关系（弱关系反倒变成了强关系），没达到上述描述的作用，
% 我将此代码进行了修改，达到了描述中的作用 ---曹磊，2023/6/9

if (nargin < 2)
    K = 0;  % 没输入K，则默认用所有存在的关系
end

N = size(CMat,1);
CAbs = abs(CMat);  % 对系数矩阵取绝对值

[~,Ind] = sort( CAbs,1,'descend' );

if (K == 0)
    for i = 1:N
        CAbs(:,i) = CAbs(:,i) ./ (CAbs(Ind(1,i),i)+eps);  % 对每列进行MinMax归一化[0,1]
    end
else
    % 修改了这一块，声明了新的变量C，起到了只取K个强关系
    C = zeros(N);
    for i = 1:N
        for j = 1:K
            C(Ind(j,i),i) = CAbs(Ind(j,i),i) ./ (CAbs(Ind(1,i),i)+eps);  % 只取最强的K个表示系数，对每列进行MinMax归一化[0,1]
        end
    end
    CAbs = C;
end

CKSym = CAbs + CAbs';  % 对称
end