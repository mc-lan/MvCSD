addpath([pwd '/dataset'])
addpath([pwd '/ClusteringMeasure'])
addpath([pwd '/code'])
clear; clc;
load('ORL_mtv.mat')
alpha = 0.7;
dim = 150;
T = 10;
for i = 1:30
    [acc(i),nmi(i),ar(i),p(i),f(i),r(i)] = MvCSD(X,gt,alpha,dim,T);
end
ACC_mean = mean(acc); ACC_std = std(acc);
NMI_mean = mean(nmi); NMI_std = std(nmi);
AR_mean = mean(ar); AR_std = std(ar);
P_mean = mean(p); P_std = std(p);
F_mean = mean(f); F_srd = std(f);
R_mean = mean(r); R_std = std(r);