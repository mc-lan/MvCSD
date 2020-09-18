addpath([pwd '/dataset'])
addpath([pwd '/code'])
clear; clc;
load('ORL_mtv.mat')
NumOfLabled = 1; % 1:10%; 2:20% ...
alpha = 0.7;
lambda = 0.1;
dim = 150;
T = 10;
for i=1:30
    [Xl,Xu,Yl,Yu] = DataCreate(X,gt,NumOfLabled);
    acc(i) = MvCSD_SS(Xl,Xu,Yl,Yu,alpha,lambda,dim,T);
end
ACC_mean = mean(acc); ACC_std = std(acc);
