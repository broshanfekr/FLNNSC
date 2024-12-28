clear;
clc;

%%
clear
warning off all
data_dir = '../dataset';
addpath(genpath(data_dir));
addpath(genpath('../util'))
addpath(genpath('../performance'))

%load USPS_1000.mat
load cifar100.mat
dataset_name = "cifar100";

%%
FEA = zeros(size(fea, 1), 32*32);
for i = 1:size(fea, 1)
   tmp_img = fea(i, :, :, :);
   tmp_img = squeeze(tmp_img);
   tmp_img = rgb2gray(tmp_img);
   tmp_img = reshape(double(tmp_img), 1, 32*32);
   FEA(i, :) = tmp_img(1,:);
end
fea = FEA;
fea = fea';
gnd = double(gnd');
X = X';

data = X;
nCluster = max(gnd)+1;
%% 
dim = nCluster * 6;  % 降维至nCluster*6；
dim = 768;
% PCA Projection
[ eigvector , eigvalue ] = PCA(data);
data = eigvector(:,1:dim)'*data;
for jj = 1 : size(data,2)
   data(:,jj) = data(:,jj)/norm(data(:,jj));  % 对data的每一列（每个样本）进行归一化
end

%% flnnsc
fprintf('flnnsc');
fprintf('Number of subjects：%d',nCluster);
para.alpha = 1e2;
para.beta = 1e1;
para.mu = 1e-4;
para.knn = 4;
para.elpson = 0.01;
para
maxIter = 10;
tol = 1e-4;
aff_type = 'J2'
gamma = 2


% tic
% [Z,W,H,obj] = FLNNSC(data,para,maxIter,tol);
% toc
% Aff = get_Aff(Z,data,aff_type,gamma);
% % Normalize each column of the affinity matrix
% Aff2 = Aff;
% for j = 1 : size(Aff,2)
%    Aff2(:,j) = Aff(:,j)/(max(abs(Aff(:,j)))+eps);    
% end
% [groups] = clu_ncut(Aff2,nCluster);  % 进行谱聚类
% [metrics.ca,metrics.nmi,metrics.ar,metrics.f1,~,~] = compute_metrics(gnd,groups);
% metrics



% %
Range = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4];
if exist(dataset_name+"_clip.mat")
    load(dataset_name+"_clip.mat")
else
    CA = zeros(length(Range));
    NMI = zeros(length(Range));
    AR = zeros(length(Range));
    F1= zeros(length(Range));
end

for i_alpha = 1:length(Range)
    para.alpha = Range(i_alpha);
    for i_beta = 1:length(Range)
        para.beta = Range(i_beta);
        if CA(i_alpha,i_beta) ~= 0
            continue
        end
        [Z,W,H,obj] = FLNNSC(data,para,maxIter,tol);
        Aff = get_Aff(Z,data,aff_type,gamma);
        % Normalize each column of the affinity matrix
        Aff2 = Aff;
        for j = 1 : size(Aff,2)
            Aff2(:,j) = Aff(:,j)/(max(abs(Aff(:,j)))+eps);    
        end
        [groups] = clu_ncut(Aff2,nCluster);  % Perform spectral clustering
        [CA(i_alpha,i_beta),NMI(i_alpha,i_beta),...
        AR(i_alpha,i_beta),F1(i_alpha,i_beta),~,~] = compute_metrics(gnd,groups);
        save(dataset_name+"_clip.mat", 'CA', "F1", "AR", "NMI")
    end
end

save(dataset_name+"_clip.mat", 'CA', "F1", "AR", "NMI")

