clear;
clc;

%%
clear
warning off all
data_dir = '../dataset';
addpath(genpath(data_dir));
addpath(genpath('../util'))
addpath(genpath('../performance'))

addpath("measure");
addpath("util");

%load USPS_1000.mat
dataset_name = "stl10";
data_path = strcat(dataset_name, ".mat");
load(data_path)

%%
FEA = zeros(size(fea, 1), size(fea, 2)*size(fea, 3));
for i = 1:size(fea, 1)
   tmp_img = fea(i, :, :, :);
   tmp_img = squeeze(tmp_img);
   tmp_img = rgb2gray(tmp_img);
   tmp_img = reshape(double(tmp_img), 1, size(fea, 2)*size(fea, 3));
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

%%
Range = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4];
if exist(dataset_name+"_clip.mat")
    load(dataset_name+"_clip.mat")

    % Find the maximum value and its linear index
    [maxValue, linearIdx] = max(NMI(:));
    % Convert linear index to row and column indices
    [row, col] = ind2sub(size(NMI), linearIdx);
    fprintf("ACC is: %.2f    NMI is: %.2f    ARI is: %.2f\n", CA(row, col)*100, NMI(row, col)*100, AR(row, col)*100)

    para.alpha = Range(row);
    para.beta = Range(col);
else
    CA = zeros(length(Range));
    NMI = zeros(length(Range));
    AR = zeros(length(Range));
    F1= zeros(length(Range));
end

%%



tic
[Z,W,H,obj] = FLNNSC(data,para,maxIter,tol);
toc
Aff = get_Aff(Z,data,aff_type,gamma);
% Normalize each column of the affinity matrix
Aff2 = Aff;
for j = 1 : size(Aff,2)
   Aff2(:,j) = Aff(:,j)/(max(abs(Aff(:,j)))+eps);    
end
[groups] = clu_ncut(Aff2,nCluster);  % 进行谱聚类
[metrics.ca,metrics.nmi,metrics.ar,metrics.f1,~,~] = compute_metrics(gnd,groups);
metrics



