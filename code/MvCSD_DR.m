function [P] = MvCSD_DR(X,alpha,dim,T,iter_out,DEBUG)
% The code is written by Mengcheng Lan,
% if you have any problems, please don't hesitate to contact me: lanmengchengds@gmail.com 
% If you find the code is useful, please cite the following reference:
% Min Meng, Mengcheng Lan, Jun Yu, Jigang Wu, 
% Multi-View Consensus Structure Discovery [J], 
% IEEE Transactions on Cybernetics, 2020.

if nargin<6
    DEBUG = 0;
end
if nargin<5
    iter_out = 3;
end
if nargin<4
    T = 30;
end

numOfView = length(X);
for i=1:numOfView
    X{i} = X{i}./(repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1));
    options.ReducedDim = dim;
    [P{i},~] = PCA1(X{i}', options);
end
for t=1:iter_out
    [Z,~] = LRS(X,P,alpha,dim,T,DEBUG);
    for i=1:numOfView
        temp = P{i}'*(X{i}-X{i}*Z);
        D = diag(0.5./sqrt(sum(temp.*temp))+eps);
        S = (X{i}-X{i}*Z)*D*(X{i}-X{i}*Z)';
        S = (S+S')/2;
        [Pall,DS]=eig(S);
        [ds,ind]=sort(diag(DS),'ascend');
        Pall=Pall(:,ind);
        [ind2]=find(ds>10^-5);
        P{i}=Pall(:,ind2(1:dim));
    end
end
end