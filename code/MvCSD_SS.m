function [ acc] = MvCSD_SS(Xl,Xu,Yl,Yu,alpha,lambda,k,T,iter_out,DEBUG)
% The code is written by Mengcheng Lan,
% if you have any problems, please don't hesitate to contact me: lanmengchengds@gmail.com 
% If you find the code is useful, please cite the following reference:
% Min Meng, Mengcheng Lan, Jun Yu, Jigang Wu, 
% Multi-View Consensus Structure Discovery [J], 
% IEEE Transactions on Cybernetics, 2020.

if nargin<10
    DEBUG = 0;
end
if nargin<9
    iter_out = 3;
end

numOfView = length(Xl);
nl = size(Yl,1);
nu = size(Yu,1);
class = length(unique(Yl));
options.ReducedDim = k;
for i=1:numOfView
    X{i} = [Xl{i} Xu{i}];
    X{i} = X{i}./(repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1));
    [P{i},~] = PCA1(X{i}', options);
end

knn_model = fitcknn(Xl{1}',Yl,'NumNeighbors',1);
Y_tar_pseudo = knn_model.predict(Xu{1}');

YYl = zeros(nl,class);
for i=1:nl
    c = Yl(i);
    YYl(i,c) = 1;
end

YYu = zeros(nu,class);
for i=1:nu
    c = Y_tar_pseudo(i);
    YYu(i,c) = 1;
end
F0 = [YYl;YYu];

for t=1:iter_out
    if t==1
        F = F0;
    end
    [Z,F,~]= LRS_SS(Xl,Xu,Yl,Yu,P,F,alpha,lambda,k,T,DEBUG);
    
    if t< iter_out
        for i=1:numOfView
            temp = P{i}'*(X{i}-X{i}*Z);
            D = diag(0.5./sqrt(sum(temp.*temp))+eps);
            S = (X{i}-X{i}*Z)*D*(X{i}-X{i}*Z)';
            S = (S+S')/2;
            [Pall,DS]=eig(S);
            [ds,ind]=sort(diag(DS),'ascend');
            Pall=Pall(:,ind);
            [ind2]=find(ds>10^-5);
            P{i}=Pall(:,ind2(1:k));
        end
    end
    Fu = F(nl+1:end,:);
end
[~,index] = max(Fu,[],2);
acc = length(find(Yu == index))/length(Yu);
end