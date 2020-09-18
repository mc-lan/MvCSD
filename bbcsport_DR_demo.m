addpath([pwd '/dataset'])
addpath([pwd '/code'])
clear; clc;
load('bbcsport_2view.mat')
alpha = 0.7;
dim = 20;
T = 30;
per = 0.4;

V = length(X);
P = MvCSD_DR(X,alpha,dim,T);
Y = [];
for i=1:V
    Y = [Y;P{i}'*X{i}];
end
for k=1:30
    [train_index,test_index] = cal_index_of_trainAndtest(gt,per);
    acc(k) = cal_acc(Y,gt,train_index,test_index);
end
acc = mean(acc);