%% 加载数据
clear
warning off all
data_dir = '../dataset';
addpath(genpath(data_dir));
addpath(genpath('../util'))
addpath(genpath('../performance'))
load Yale.mat
%% 预处理---PCA降维
nCluster = 15;
dim = nCluster * 6;  % 降维至nCluster*6；
%% PCA Projection
[ eigvector , eigvalue ] = PCA(fea) ;
data = eigvector(:,1:dim)'*fea;
for jj = 1 : size(data,2)
   data(:,jj) = data(:,jj)/norm(data(:,jj));  % 对data的每一列（每个样本）进行归一化
end
%% flnnsc单次实验
para.alpha = 1e2;
para.beta = 1e0;
para.mu = 1e-4;
para.knn = 4;
para.elpson = 0.001;
maxIter = 10;
tol = 1e-2;
aff_type = 'J2';
gamma = 2;

tic
[Z,W,H,obj] = FLNNSC(data,para,maxIter,tol);
toc
Aff = get_Aff(Z,data,aff_type,gamma);
% 对亲和矩阵的每列进行归一化
Aff2 = Aff;
for j = 1 : size(Aff,2)
    Aff2(:,j) = Aff(:,j)/(max(abs(Aff(:,j)))+eps);    
end
[groups] = clu_ncut(Aff2,nCluster);  % 进行谱聚类
[metrics.ca,metrics.nmi,metrics.ar,metrics.f1,~,~] = compute_metrics(gnd,groups);
metrics